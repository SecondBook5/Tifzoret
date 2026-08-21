localrules: all

rule all:
    input:
        str(RESULTS / "manifest.json")
    default_target: True

rule materialize_inputs:
    input:
        sources=SOURCE_FILES,
        samples=SOURCE_SAMPLES,
        contrasts=SOURCE_CONTRASTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "materialize_inputs.py")
    output:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        manifest=INPUT_MANIFEST,
        **COMPANION_OUTPUTS,
        **ABUNDANCE_OUTPUTS
    log:
        str(INPUTS / "materialize.log")
    params:
        companion=COMPANION_ARGS,
        abundance=ABUNDANCE_ARGS
    threads:
        CONFIG.get("counting", {}).get("threads", 1)
    conda:
        CORE_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --counts {output.counts:q} "
        "--samples {output.samples:q} --annotation {output.annotation:q} "
        "--contrasts {output.contrasts:q} {params.companion} {params.abundance} "
        "--manifest {output.manifest:q} --threads {threads} > {log:q} 2>&1"

rule study_qc:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "qc.R"),
        utils=UTILS_R
    output:
        vst=qc("objects/vst.rds"),
        expression=qc("tables/vst_expression.tsv"),
        pca_coordinates=qc("tables/pca_coordinates.tsv"),
        pca_variance=qc("tables/pca_variance.tsv"),
        correlation=qc("tables/sample_correlation.tsv"),
        library_metrics=qc("tables/library_metrics.tsv"),
        density_table=qc("tables/expression_density_displayed.tsv"),
        distance_table=qc("tables/sample_distance.tsv"),
        variable_table=qc("tables/variable_gene_heatmap_displayed.tsv"),
        outlier_table=qc("tables/outlier_distances.tsv"),
        pca_correlation_layout=qc("tables/pca_correlation_layout.json"),
        pca_pdf=qc("figures/pca.pdf"),
        pca_png=qc("figures/pca.png"),
        correlation_pdf=qc("figures/sample_correlation.pdf"),
        correlation_png=qc("figures/sample_correlation.png"),
        pca_correlation_pdf=qc("figures/pca_correlation.pdf"),
        pca_correlation_png=qc("figures/pca_correlation.png"),
        metrics_pdf=qc("figures/library_metrics.pdf"),
        metrics_png=qc("figures/library_metrics.png"),
        density_pdf=qc("figures/expression_density.pdf"),
        density_png=qc("figures/expression_density.png"),
        distance_pdf=qc("figures/sample_distance.pdf"),
        distance_png=qc("figures/sample_distance.png"),
        variable_pdf=qc("figures/variable_gene_heatmap.pdf"),
        variable_png=qc("figures/variable_gene_heatmap.png"),
        outlier_pdf=qc("figures/sample_outliers.pdf"),
        outlier_png=qc("figures/sample_outliers.png"),
        overview_pdf=qc("figures/qc_overview.pdf"),
        overview_png=qc("figures/qc_overview.png"),
        summary=qc("qc_summary.json")
    log:
        qc("logs/qc.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--counts {input.counts:q} --samples {input.samples:q} --annotation {input.annotation:q} "
        "--outdir {RESULTS}/qc > {log:q} 2>&1"

rule study_batch:
    input:
        vst=qc("objects/vst.rds"),
        samples=SAMPLES,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "batch.R"),
        utils=UTILS_R
    output:
        corrected=batch("tables/batch_corrected_expression.tsv"),
        pca_coordinates=batch("tables/batch_pca_coordinates.tsv"),
        distance=batch("tables/batch_sample_distance.tsv"),
        pca_pdf=batch("figures/batch_pca.pdf"),
        pca_png=batch("figures/batch_pca.png"),
        distance_pdf=batch("figures/batch_sample_distance.pdf"),
        distance_png=batch("figures/batch_sample_distance.png"),
        summary=batch("batch_summary.json")
    log:
        batch("logs/batch.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--vst {input.vst:q} --samples {input.samples:q} "
        "--outdir {RESULTS}/batch > {log:q} 2>&1"

rule study_deconvolution:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        signature=str(PROJECT.deconvolution_signature) if PROJECT.deconvolution_signature else [],
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "deconvolution.R"),
        utils=UTILS_R
    output:
        fractions=deconvolution("tables/cell_fractions.tsv"),
        fractions_pdf=deconvolution("figures/cell_fractions.pdf"),
        fractions_png=deconvolution("figures/cell_fractions.png"),
        summary=deconvolution("deconvolution_summary.json")
    log:
        deconvolution("logs/deconvolution.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--counts {input.counts:q} --samples {input.samples:q} --annotation {input.annotation:q} "
        "--signature {input.signature:q} --outdir {RESULTS}/deconvolution > {log:q} 2>&1"

rule study_consensus:
    input:
        de=expand(analysis("de", "tables/de_results.tsv"), contrast_id=PAIRWISE_CONTRAST_IDS),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "consensus.py")
    output:
        membership=comparison("tables/consensus_membership.tsv"),
        genes=comparison("tables/consensus_genes.tsv"),
        overlap=comparison("tables/contrast_overlap.tsv"),
        intersections=comparison("tables/consensus_intersections_displayed.tsv"),
        upset_pdf=comparison("figures/consensus_upset.pdf"),
        upset_png=comparison("figures/consensus_upset.png"),
        heatmap_pdf=comparison("figures/consensus_sign_heatmap.pdf"),
        heatmap_png=comparison("figures/consensus_sign_heatmap.png"),
        summary=comparison("consensus_summary.json")
    log:
        comparison("logs/consensus.log")
    conda:
        NETWORK_ENV
    shell:
        "python {input.script} --project-config {input.config:q} --de {input.de:q} "
        "--outdir {RESULTS}/comparison > {log:q} 2>&1"

rule study_variance_partition:
    input:
        vst=qc("objects/vst.rds"),
        samples=SAMPLES,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "variancepartition.R"),
        utils=UTILS_R
    output:
        fractions=variance_partition("tables/variance_fractions.tsv"),
        summary_table=variance_partition("tables/variance_summary.tsv"),
        pdf=variance_partition("figures/variance_partition.pdf"),
        png=variance_partition("figures/variance_partition.png"),
        summary=variance_partition("variance_partition_summary.json")
    log:
        variance_partition("logs/variance_partition.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--vst {input.vst:q} --samples {input.samples:q} "
        "--outdir {RESULTS}/variance_partition > {log:q} 2>&1"

rule contrast_de:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "de.R"),
        utils=UTILS_R
    output:
        dds=analysis("de", "objects/deseq2.rds"),
        results=analysis("de", "tables/de_results.tsv"),
        volcano_table=analysis("de", "tables/volcano_displayed.tsv"),
        heatmap_table=analysis("de", "tables/de_heatmap_displayed.tsv"),
        ma_table=analysis("de", "tables/ma_displayed.tsv"),
        pvalue_table=analysis("de", "tables/pvalue_distribution_displayed.tsv"),
        lfc_table=analysis("de", "tables/lfc_distribution_displayed.tsv"),
        de_pca_table=analysis("de", "tables/de_pca_coordinates.tsv"),
        de_pca_ellipses=analysis("de", "tables/de_pca_ellipses.tsv"),
        volcano_pdf=analysis("de", "figures/volcano.pdf"),
        volcano_png=analysis("de", "figures/volcano.png"),
        heatmap_pdf=analysis("de", "figures/de_heatmap.pdf"),
        heatmap_png=analysis("de", "figures/de_heatmap.png"),
        ma_pdf=analysis("de", "figures/ma.pdf"),
        ma_png=analysis("de", "figures/ma.png"),
        pvalue_pdf=analysis("de", "figures/pvalue_distribution.pdf"),
        pvalue_png=analysis("de", "figures/pvalue_distribution.png"),
        lfc_pdf=analysis("de", "figures/lfc_distribution.pdf"),
        lfc_png=analysis("de", "figures/lfc_distribution.png"),
        de_pca_pdf=analysis("de", "figures/de_pca.pdf"),
        de_pca_png=analysis("de", "figures/de_pca.png"),
        overview_pdf=analysis("de", "figures/de_overview.pdf"),
        overview_png=analysis("de", "figures/de_overview.png"),
        summary=analysis("de", "de_summary.json")
    log:
        analysis("de", "logs/de.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--counts {input.counts:q} --samples {input.samples:q} --annotation {input.annotation:q} "
        "--contrasts {input.contrasts:q} "
        "--contrast-id {wildcards.contrast_id:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/de > {log:q} 2>&1"

rule contrast_omnibus:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "omnibus.R"),
        utils=UTILS_R
    output:
        dds=analysis("omnibus", "objects/deseq2_lrt.rds"),
        results=analysis("omnibus", "tables/omnibus_results.tsv"),
        pvalue_table=analysis("omnibus", "tables/pvalue_distribution_displayed.tsv"),
        heatmap_table=analysis("omnibus", "tables/omnibus_heatmap_displayed.tsv"),
        pvalue_pdf=analysis("omnibus", "figures/pvalue_distribution.pdf"),
        pvalue_png=analysis("omnibus", "figures/pvalue_distribution.png"),
        heatmap_pdf=analysis("omnibus", "figures/omnibus_heatmap.pdf"),
        heatmap_png=analysis("omnibus", "figures/omnibus_heatmap.png"),
        summary=analysis("omnibus", "omnibus_summary.json")
    log:
        analysis("omnibus", "logs/omnibus.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--counts {input.counts:q} --samples {input.samples:q} --annotation {input.annotation:q} "
        "--contrasts {input.contrasts:q} "
        "--contrast-id {wildcards.contrast_id:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/omnibus > {log:q} 2>&1"
