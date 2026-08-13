rule contrast_publication:
    input:
        vst=qc("objects/vst.rds"),
        de=analysis("de", "tables/de_results.tsv"),
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        panels=str(PROJECT.hypothesis_panels),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "publication.R"),
        utils=UTILS_R
    output:
        PUBLICATION_PATTERNS
    log:
        analysis("publication", "logs/publication.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} --samples {input.samples:q} "
        "--annotation {input.annotation:q} --contrasts {input.contrasts:q} --contrast-id {wildcards.contrast_id:q} "
        "--vst {input.vst:q} --de {input.de:q} --panels {input.panels:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/publication > {log:q} 2>&1"

rule assemble_figure:
    input:
        panels=figure_panel_inputs,
        config=str(CONFIG_PATH),
        recipe=str(PROJECT.figure_recipe) if PROJECT.figure_recipe else []
    output:
        pdf=str(RESULTS / "publication" / "{figure_set}" / "assembled" / "{figure_set}.pdf"),
        png=str(RESULTS / "publication" / "{figure_set}" / "assembled" / "{figure_set}.png"),
        metadata=str(RESULTS / "publication" / "{figure_set}" / "assembled" / "assembly.json"),
        panel_index=str(RESULTS / "publication" / "{figure_set}" / "panels" / "index.json")
    log:
        str(RESULTS / "publication" / "{figure_set}" / "assembled" / "assembly.log")
    conda:
        CORE_ENV
    shell:
        "python {WORKFLOW_ROOT}/scripts/assemble.py --project-config {input.config:q} "
        "--results {RESULTS:q} --figure-set {wildcards.figure_set:q} "
        "--pdf {output.pdf:q} --png {output.png:q} --metadata {output.metadata:q} "
        "--panel-index {output.panel_index:q} > {log:q} 2>&1"

rule assemble_publication:
    input:
        ASSEMBLY_OUTPUTS

rule front_door_artifacts:
    input:
        artifacts=ANALYSIS_OUTPUTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "front_door.py")
    output:
        figures=str(RESULTS / "figures" / "index.json"),
        tables=str(RESULTS / "tables" / "index.json")
    log:
        str(RESULTS / ".cache" / "logs" / "front_door.log")
    conda:
        CORE_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --results {RESULTS:q} "
        "--figures-index {output.figures:q} --tables-index {output.tables:q} > {log:q} 2>&1"
