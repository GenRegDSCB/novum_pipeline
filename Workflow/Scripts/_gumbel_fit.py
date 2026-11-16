"""Gumbel-fit quality control and the per-miRNA theta(L) power law.

Stdlib only: `rule rnacalibrate` runs in `rnahybrid.yaml`, which has no pandas or numpy.

The model:  theta(L) = C * [ln(m*n)] ** -alpha   (m = target length, n = query length),
a straight line in logs, so alpha is closed-form weighted least squares over a miRNA's own
anchors. Weights are 1/SE^2 with SE(ln theta) = 1.053/sqrt(N), from the measured
SE(theta_hat) = 1.053*theta/sqrt(N).

Tier 1 is not heuristic: every condition is a literal branch in RNAcalibrate's
`numerical.c`.
"""

import math
from collections import namedtuple

# The asymptotic Gumbel-scale coefficient, confirmed across repeated draws at identical args.
SE_COEFFICIENT = 1.053

# Two points fix a line exactly, leaving no residual to test.
MIN_POINTS = 3

FitPoint = namedtuple("FitPoint", "target_length query_length theta sample_size")
AlphaFit = namedtuple("AlphaFit", "alpha log_c kept dropped residual_spread")


#----- Tier 3 has no usable curve. Its own class so the caller falls back by TYPE, not by
#      matching a message -----#
class UnusableCurveError(RuntimeError):
    pass


#----- SE of ln(theta) for one fit -----#
def se_log_theta(sample_size):
    return SE_COEFFICIENT / math.sqrt(sample_size)


#----- Tier 1 hard rejects, each a literal branch of the RNAcalibrate source. None = sound.
#      ORDER MATTERS: non-finite precedes the sign checks, since NaN <= 0 is False -----#
def classify_fit(sample_size, xi, theta):
    # Column 2 is `end-start+1`, a difference of bin indices (can go negative), and the
    # zeros branch is `if (end<=start)` -- i.e. N <= 1, not N < 0.
    if sample_size <= 1:
        return "degenerate_sample"
    if not (math.isfinite(xi) and math.isfinite(theta)):
        return "non_finite"
    # The dangerous class: finite, passes any isnan guard. Measured: xi=13970/theta=-2022,
    # xi=-2070/theta=300, xi=-3.33/theta=0.803.
    if theta <= 0:
        return "non_positive_theta"
    if xi <= 0:
        return "non_positive_xi"
    return None


#----- theta off a fitted curve at any target length -----#
def predict_theta(alpha, log_c, target_length, query_length):
    return math.exp(log_c - alpha * math.log(math.log(target_length * query_length)))


#----- Observed theta over predicted; 1.0 is exactly on the curve -----#
def residual_ratio(point, alpha, log_c):
    return point.theta / predict_theta(alpha, log_c, point.target_length, point.query_length)


#----- One weighted least-squares pass; fit_alpha_robust layers rejection on top -----#
def fit_alpha(points):
    if len(points) < MIN_POINTS:
        raise UnusableCurveError(
            f"Fitting alpha needs at least {MIN_POINTS} anchors, got {len(points)}."
        )

    sum_w = sum_x = sum_y = sum_xx = sum_xy = 0.0
    for point in points:
        x = math.log(math.log(point.target_length * point.query_length))
        y = math.log(point.theta)
        weight = point.sample_size / SE_COEFFICIENT ** 2
        sum_w += weight
        sum_x += weight * x
        sum_y += weight * y
        sum_xx += weight * x * x
        sum_xy += weight * x * y

    denominator = sum_w * sum_xx - sum_x * sum_x
    if denominator == 0:
        raise UnusableCurveError(
            "Fitting alpha needs anchors at more than one length; every point shares "
            "the same ln(ln(m*n))."
        )

    slope = (sum_w * sum_xy - sum_x * sum_y) / denominator
    log_c = (sum_y * sum_xx - sum_x * sum_xy) / denominator
    alpha = -slope

    ratios = [residual_ratio(point, alpha, log_c) for point in points]
    spread = max(ratios) / min(ratios) - 1.0
    return AlphaFit(alpha=alpha, log_c=log_c, kept=list(points), dropped=[],
                    residual_spread=spread)


#----- Tier 3: fit, drop beyond n_sigma of the KNOWN per-fit SE, refit. One pass, not an
#      iteration -- the threshold comes from N, not the residual scale, so an outlier cannot
#      mask itself. Raises rather than return a nonsense curve -----#
def fit_alpha_robust(points, n_sigma=3.0):
    first = fit_alpha(points)

    kept, dropped = [], []
    for point in points:
        z = math.log(residual_ratio(point, first.alpha, first.log_c)) / se_log_theta(point.sample_size)
        (dropped if abs(z) > n_sigma else kept).append(point)

    if dropped:
        if len(kept) < MIN_POINTS:
            raise UnusableCurveError(
                f"Only {len(kept)} anchor(s) survived {n_sigma}-sigma rejection; alpha "
                f"needs at least {MIN_POINTS}."
            )
        refit = fit_alpha(kept)
        result = refit._replace(dropped=dropped)
    else:
        result = first

    # theta must FALL with length; a flat or rising curve is not describing the null.
    if result.alpha <= 0:
        raise UnusableCurveError(
            f"Fitted alpha is {result.alpha:.4f}; theta must decrease with target length, "
            "so this curve cannot be extrapolated from."
        )
    return result


#----- Value at which cumulative weight first reaches half. Used for xi: it drifts only a few
#      percent across a ladder but is noisy per fit, so the survivors' median wins -----#
def weighted_median(values, weights):
    if not values:
        raise ValueError("weighted_median needs at least one value.")
    ordered = sorted(zip(values, weights))
    half = sum(weights) / 2.0
    running = 0.0
    for value, weight in ordered:
        running += weight
        if running >= half:
            return value
    # Unreachable on non-negative weights, which is all a sample_size can be; here so a
    # negative one degrades to the largest value instead of returning None.
    return ordered[-1][0]
