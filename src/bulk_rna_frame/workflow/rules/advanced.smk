rule contrast_sva:
    input: dds=analysis("de", "objects/deseq2.rds"), vst=qc("objects/vst.rds"), contrasts=CONTRASTS, config=str(CONFIG_PATH), script=str(WORKFLOW_ROOT / "scripts" / "sva.R"), utils=UTILS_R
    output: SVA_PATTERNS
    log: analysis("advanced/sva", "logs/sva.log")
    conda: R_ENV
    shell: "Rscript --vanilla {input.script} --project-config {input.config:q} --dds {input.dds:q} --vst {input.vst:q} --contrasts {input.contrasts:q} --contrast-id {wildcards.contrast_id:q} --outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/advanced/sva > {log:q} 2>&1"

rule contrast_wgcna:
    input: vst=qc("objects/vst.rds"), samples=SAMPLES, annotation=ANNOTATION, contrasts=CONTRASTS, config=str(CONFIG_PATH), script=str(WORKFLOW_ROOT / "scripts" / "wgcna.R"), utils=UTILS_R
    output: WGCNA_PATTERNS
    log: analysis("advanced/wgcna", "logs/wgcna.log")
    conda: R_ENV
    shell: "Rscript --vanilla {input.script} --project-config {input.config:q} --vst {input.vst:q} --samples {input.samples:q} --annotation {input.annotation:q} --contrasts {input.contrasts:q} --contrast-id {wildcards.contrast_id:q} --outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/advanced/wgcna > {log:q} 2>&1"

rule contrast_mediation:
    input: scores=analysis("pathways", "tables/gsva_scores.tsv"), samples=SAMPLES, contrasts=CONTRASTS, config=str(CONFIG_PATH), script=str(WORKFLOW_ROOT / "scripts" / "mediation.R"), utils=UTILS_R
    output: inputs=analysis("advanced/mediation", "tables/mediation_inputs.tsv"), results=analysis("advanced/mediation", "tables/mediation_results.tsv"), summary=analysis("advanced/mediation", "mediation_summary.json")
    log: analysis("advanced/mediation", "logs/mediation.log")
    conda: R_ENV
    shell: "Rscript --vanilla {input.script} --project-config {input.config:q} --scores {input.scores:q} --samples {input.samples:q} --contrasts {input.contrasts:q} --contrast-id {wildcards.contrast_id:q} --outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/advanced/mediation > {log:q} 2>&1"

rule contrast_mediation_power:
    input: mediation=analysis("advanced/mediation", "tables/mediation_inputs.tsv"), script=str(WORKFLOW_ROOT / "scripts" / "power.py")
    output: table=analysis("advanced/mediation", "tables/mediation_power.tsv"), summary=analysis("advanced/mediation", "mediation_power_summary.json")
    log: analysis("advanced/mediation", "logs/power.log")
    conda: CORE_ENV
    shell: "python {input.script} --inputs {input.mediation:q} --outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/advanced/mediation > {log:q} 2>&1"

rule contrast_multilayer:
    input: grn=analysis("regulators", "tables/grn_edges_displayed.tsv"), hubs=analysis("advanced/wgcna", "tables/wgcna_hubs.tsv"), string_up=analysis("networks", "tables/string_up_edges.tsv"), string_down=analysis("networks", "tables/string_down_edges.tsv"), script=str(WORKFLOW_ROOT / "scripts" / "multilayer.py")
    output: MULTILAYER_PATTERNS
    log: analysis("advanced/multilayer", "logs/multilayer.log")
    conda: NETWORK_ENV
    shell: "python {input.script} --grn-edges {input.grn:q} --wgcna-hubs {input.hubs:q} --string-up {input.string_up:q} --string-down {input.string_down:q} --contrast-id {wildcards.contrast_id:q} --outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/advanced/multilayer > {log:q} 2>&1"
