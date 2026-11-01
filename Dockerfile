FROM mambaorg/micromamba:2.3.2

USER root
# micromamba writes a temporary pip requirements file into the working directory, so the unprivileged user must own it.
RUN mkdir -p /opt/tifzoret && chown $MAMBA_USER:$MAMBA_USER /opt/tifzoret
WORKDIR /opt/tifzoret
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/tifzoret
USER $MAMBA_USER
RUN micromamba install --yes --name base --file environment.yaml \
    && micromamba clean --all --yes

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "tifzoret"]
CMD ["--help"]

