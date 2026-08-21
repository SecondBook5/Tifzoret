rule report_html:
    input:
        artifacts=REPORT_INPUTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "report.py")
    output:
        REPORT
    log:
        str(RESULTS / ".cache" / "logs" / "report.log")
    conda:
        CORE_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --results {RESULTS:q} "
        "--output {output:q} > {log:q} 2>&1"

rule release_manifest:
    input:
        artifacts=FINAL_OUTPUTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "manifest.py")
    output:
        str(RESULTS / "manifest.json")
    log:
        str(RESULTS / "manifest.log")
    params:
        results=str(RESULTS)
    conda:
        CORE_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --results {params.results:q} --output {output:q} > {log:q} 2>&1"
