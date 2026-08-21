FROM mambaorg/micromamba:2.3.2

USER root
WORKDIR /opt/tifzoret
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/tifzoret
USER $MAMBA_USER
RUN micromamba install --yes --name base --file environment.yaml \
    && micromamba clean --all --yes

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "tifzoret"]
CMD ["--help"]

