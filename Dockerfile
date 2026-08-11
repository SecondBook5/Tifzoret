FROM mambaorg/micromamba:2.3.2

USER root
WORKDIR /opt/bulk-rna-frame
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/bulk-rna-frame
USER $MAMBA_USER
RUN micromamba install --yes --name base --file environment.yaml \
    && micromamba clean --all --yes

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "bulk-rna"]
CMD ["--help"]

