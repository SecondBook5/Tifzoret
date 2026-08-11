rule contrast_pathways:
    input:
        de=analysis("de", "tables/de_results.tsv"),
        vst=qc("objects/vst.rds"),
        samples=SAMPLES,
        annotation=ANNOTATION,
        gmt=GMT,
        resources=RESOURCE_TABLE,
        panels=str(PROJECT.hypothesis_panels) if PROJECT.hypothesis_panels else [],
        config=str(CONFIG_PATH),
        contrasts=CONTRASTS,
        script=str(WORKFLOW_ROOT / "scripts" / "pathways.R"),
        utils=UTILS_R
    output:
        fgsea=analysis("pathways", "tables/fgsea.tsv"),
        ora=analysis("pathways", "tables/ora.tsv"),
        ora_displayed=analysis("pathways", "tables/ora_displayed.tsv"),
        gsva_scores=analysis("pathways", "tables/gsva_scores.tsv"),
        gsva_differential=analysis("pathways", "tables/gsva_differential.tsv"),
        gsva_displayed=analysis("pathways", "tables/gsva_heatmap_displayed.tsv"),
        gsea_displayed=analysis("pathways", "tables/gsea_curves_displayed.tsv"),
        ora_pdf=analysis("pathways", "figures/ora_bidirectional.pdf"),
        ora_png=analysis("pathways", "figures/ora_bidirectional.png"),
        gsva_pdf=analysis("pathways", "figures/gsva_heatmap.pdf"),
        gsva_png=analysis("pathways", "figures/gsva_heatmap.png"),
        gsea_pdf=analysis("pathways", "figures/gsea_curves.pdf"),
        gsea_png=analysis("pathways", "figures/gsea_curves.png"),
        summary=analysis("pathways", "pathways_summary.json")
    log:
        analysis("pathways", "logs/pathways.log")
    params:
        panels_option=("--panels " + shlex.quote(str(PROJECT.hypothesis_panels))) if PROJECT.hypothesis_panels else ""
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--samples {input.samples:q} --annotation {input.annotation:q} "
        "--contrasts {input.contrasts:q} --gmt {input.gmt:q} --resource-table {input.resources:q} {params.panels_option} "
        "--contrast-id {wildcards.contrast_id:q} --vst {input.vst:q} --de {input.de:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/pathways > {log:q} 2>&1"

rule contrast_ontology:
    input:
        ora=analysis("pathways", "tables/ora.tsv"),
        resources=RESOURCE_TABLE,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "ontology.R"),
        utils=UTILS_R
    output:
        table=analysis("ontology", "tables/ontology.tsv"),
        displayed=analysis("ontology", "tables/ontology_displayed.tsv"),
        pdf=analysis("ontology", "figures/ontology_bidirectional.pdf"),
        png=analysis("ontology", "figures/ontology_bidirectional.png"),
        summary=analysis("ontology", "ontology_summary.json")
    log:
        analysis("ontology", "logs/ontology.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} --ora {input.ora:q} "
        "--resource-table {input.resources:q} --contrast-id {wildcards.contrast_id:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/ontology > {log:q} 2>&1"

rule contrast_composition:
    input:
        vst=qc("objects/vst.rds"),
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        signatures=str(PROJECT.cell_state_signatures) if PROJECT.cell_state_signatures else [],
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "composition.R"),
        utils=UTILS_R
    output:
        scores=analysis("composition", "tables/cell_state_scores.tsv"),
        differential=analysis("composition", "tables/cell_state_differential.tsv"),
        displayed=analysis("composition", "tables/cell_state_displayed.tsv"),
        pdf=analysis("composition", "figures/cell_state_signatures.pdf"),
        png=analysis("composition", "figures/cell_state_signatures.png"),
        summary=analysis("composition", "composition_summary.json")
    log:
        analysis("composition", "logs/composition.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} --samples {input.samples:q} "
        "--annotation {input.annotation:q} --contrasts {input.contrasts:q} --contrast-id {wildcards.contrast_id:q} "
        "--vst {input.vst:q} --signatures {input.signatures:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/composition > {log:q} 2>&1"

rule contrast_regulators:
    input:
        vst=qc("objects/vst.rds"),
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        regulon=str(PROJECT.regulon_edges) if PROJECT.regulon_edges else [],
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "regulators.R"),
        utils=UTILS_R
    output:
        edges=analysis("regulators", "tables/regulon_edges.tsv"),
        signed=analysis("regulators", "tables/dorothea_activity_scores.tsv"),
        unsigned=analysis("regulators", "tables/regulator_target_program_scores.tsv"),
        differential=analysis("regulators", "tables/regulator_differential.tsv"),
        displayed=analysis("regulators", "tables/regulator_activity_displayed.tsv"),
        pdf=analysis("regulators", "figures/regulator_activity.pdf"),
        png=analysis("regulators", "figures/regulator_activity.png"),
        summary=analysis("regulators", "regulators_summary.json")
    log:
        analysis("regulators", "logs/regulators.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} --samples {input.samples:q} "
        "--annotation {input.annotation:q} --contrasts {input.contrasts:q} --contrast-id {wildcards.contrast_id:q} "
        "--vst {input.vst:q} --outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/regulators > {log:q} 2>&1"

rule contrast_networks:
    input:
        de=analysis("de", "tables/de_results.tsv"),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "networks.py")
    output:
        NETWORK_PATTERNS
    log:
        analysis("networks", "logs/networks.log")
    conda:
        NETWORK_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --de {input.de:q} "
        "--contrast-id {wildcards.contrast_id:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/networks "
        "--cache-dir {RESULTS}/.cache/resources > {log:q} 2>&1"

rule contrast_grn:
    input:
        edges=analysis("regulators", "tables/regulon_edges.tsv"),
        regulators=analysis("regulators", "tables/regulator_differential.tsv"),
        de=analysis("de", "tables/de_results.tsv"),
        panels=str(PROJECT.hypothesis_panels) if PROJECT.hypothesis_panels else [],
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "grn.py")
    output:
        nodes=analysis("regulators", "tables/grn_nodes_displayed.tsv"),
        edges=analysis("regulators", "tables/grn_edges_displayed.tsv"),
        separation=analysis("regulators", "tables/grn_program_separation_test.tsv"),
        sectors=analysis("regulators", "tables/grn_sector_summary.tsv"),
        rectangular_pdf=analysis("regulators", "figures/grn_rectangular.pdf"),
        rectangular_png=analysis("regulators", "figures/grn_rectangular.png"),
        radial_pdf=analysis("regulators", "figures/grn_radial.pdf"),
        radial_png=analysis("regulators", "figures/grn_radial.png"),
        summary=analysis("regulators", "grn_summary.json")
    log:
        analysis("regulators", "logs/grn.log")
    conda:
        NETWORK_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --edges {input.edges:q} "
        "--regulators {input.regulators:q} --de {input.de:q} --contrast-id {wildcards.contrast_id:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/regulators > {log:q} 2>&1"

rule contrast_hypotheses:
    input:
        de=analysis("de", "tables/de_results.tsv"),
        fgsea=analysis("pathways", "tables/fgsea.tsv"),
        gsva=analysis("pathways", "tables/gsva_differential.tsv"),
        regulators=analysis("regulators", "tables/regulator_differential.tsv") if MODULES["regulators"] else [],
        claims=str(PROJECT.hypotheses),
        panels=str(PROJECT.hypothesis_panels),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "hypotheses.py")
    output:
        evidence=analysis("hypotheses", "tables/hypothesis_evidence.tsv"),
        summary_table=analysis("hypotheses", "tables/hypothesis_summary.tsv"),
        summary=analysis("hypotheses", "hypotheses_summary.json"),
        report=analysis("hypotheses", "hypotheses_report.html")
    params:
        regulator_option=lambda wildcards: (
            "--regulators " + shlex.quote(
                analysis("regulators", "tables/regulator_differential.tsv").format(
                    contrast_id=wildcards.contrast_id
                )
            ) if MODULES["regulators"] else ""
        )
    log:
        analysis("hypotheses", "logs/hypotheses.log")
    conda:
        CORE_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --contrast-id {wildcards.contrast_id:q} "
        "--de {input.de:q} --fgsea {input.fgsea:q} --gsva {input.gsva:q} {params.regulator_option} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/hypotheses > {log:q} 2>&1"
