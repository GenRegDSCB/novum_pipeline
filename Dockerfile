# syntax=docker/dockerfile:1
#----- Novum Pipeline: Snakemake plus every rule env, pre-built. linux/amd64 only: RNAhybrid and -----#
#----- libfaketime have no ARM Linux build, so Apple Silicon runs this image through Rosetta -----#
FROM --platform=linux/amd64 condaforge/miniforge3:26.9.4-0@sha256:63d72eadf0ab73ee63bad8f5cfdfe1754530415b4aeba99e466e9a6f77ceec95 AS base
RUN conda config --system --set channel_priority strict

#----- lock / locks: regenerate the linux-64 pin files from the YAMLs alone. Run with -----#
#----- docker build --platform linux/amd64 --target locks --output type=local,dest=Workflow/Envs . -----#
FROM base AS lock
COPY Workflow/Envs/*.yaml /envs/
RUN mkdir /locks && for f in /envs/*.yaml; do \
        stem=$(basename "$f" .yaml); \
        conda env create -q -f "$f" -p "/tmp/envs/$stem" && \
        conda list --explicit --md5 -p "/tmp/envs/$stem" > "/locks/$stem.linux-64.pin.txt" || exit 1; \
    done && rm -rf /tmp/envs && conda clean -afy

FROM scratch AS locks
COPY --from=lock /locks/ /

#----- runtime: what gets published. Rule envs live at a fixed prefix: Snakemake hashes it -----#
FROM base AS runtime
RUN --mount=type=bind,source=Workflow/Envs/environment.linux-64.pin.txt,target=/tmp/driver.pin.txt \
    conda create -y -p /opt/snakemake --file /tmp/driver.pin.txt && conda clean -afy

WORKDIR /opt/novum_pipeline
COPY Snakefile ./
COPY Config/config.yaml Config/
COPY Workflow/Envs/ Workflow/Envs/
COPY Workflow/Scripts/*.py Workflow/Scripts/build_plots.R Workflow/Scripts/
COPY data-prep/*.py data-prep/
COPY docker/ docker/

ENV PATH=/opt/snakemake/bin:$PATH
RUN docker/create-envs.sh && conda clean -afy && rm -rf /root/.cache

#----- Any --user works: HOME is /tmp under Docker (Apptainer keeps the host's), nothing else is -----#
#----- written at run time. No user site-packages, so a host ~/.local cannot shadow the envs -----#
ENV HOME=/tmp PYTHONNOUSERSITE=1
WORKDIR /work
ENTRYPOINT ["/opt/novum_pipeline/docker/entrypoint.sh"]

ARG VERSION=dev
ARG REVISION=unknown
LABEL org.opencontainers.image.title="novum_pipeline" \
      org.opencontainers.image.source="https://github.com/GenRegDSCB/novum_pipeline" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${REVISION}"

#----- test: runtime plus pytest and the suites, rule-env binaries appended to PATH. Never published -----#
FROM runtime AS test
RUN conda install -y -p /opt/snakemake pytest=9.1.1 && conda clean -afy
COPY data-prep/tests/ /opt/novum_pipeline/data-prep/tests/
COPY Workflow/Scripts/tests/ /opt/novum_pipeline/Workflow/Scripts/tests/
ENTRYPOINT ["bash", "-c", "PATH=\"$PATH$(printf ':%s' /opt/conda-envs/*/bin)\" exec pytest -p no:cacheprovider /opt/novum_pipeline/data-prep/tests /opt/novum_pipeline/Workflow/Scripts/tests \"$@\"", "run-tests"]

#----- release: last, so a plain `docker build .` gives the published image -----#
FROM runtime AS release
