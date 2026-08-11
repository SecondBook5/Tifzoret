rule resolve_resources:
    input:
        custom_gmt=CUSTOM_GMT,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "resources.R"),
        utils=UTILS_R
    output:
        gmt=GMT,
        table=RESOURCE_TABLE,
        receipt=RESOURCE_RECEIPT
    log:
        str(RESULTS / ".cache" / "logs" / "resources.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--custom-gmt {input.custom_gmt:q} --gmt {output.gmt:q} --table {output.table:q} "
        "--receipt {output.receipt:q} > {log:q} 2>&1"
