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
        script=str(WORKFLOW_ROOT / "scripts" / "01_inputs" / "materialize_inputs.py")
    output:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        contrasts=CONTRASTS,
        families=FAMILIES,
        estimands=ESTIMANDS,
        cell_weights=CELL_WEIGHTS,
        term_tests=TERM_TESTS,
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
        "--contrasts {output.contrasts:q} "
        "--families {output.families:q} --estimands {output.estimands:q} "
        "--cell-weights {output.cell_weights:q} --term-tests {output.term_tests:q} "
        "{params.companion} {params.abundance} "
        "--manifest {output.manifest:q} --threads {threads} > {log:q} 2>&1"

rule study_qc:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "02_qc" / "qc.R"),
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
        detected_table=qc("tables/detected_genes_displayed.tsv"),
        detected_pdf=qc("figures/detected_genes.pdf"),
        detected_png=qc("figures/detected_genes.png"),
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
        script=str(WORKFLOW_ROOT / "scripts" / "02_qc" / "batch.R"),
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
        script=str(WORKFLOW_ROOT / "scripts" / "05_composition" / "deconvolution.R"),
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
        script=str(WORKFLOW_ROOT / "scripts" / "09_synthesis" / "consensus.py")
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
        script=str(WORKFLOW_ROOT / "scripts" / "02_qc" / "variancepartition.R"),
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

rule study_factorial:
    input:
        vst_expression=qc("tables/vst_expression.tsv"),
        samples=SAMPLES,
        de_x=analysis("de", "tables/de_results.tsv").format(contrast_id=FACTORIAL_EFFECT_X),
        de_y=analysis("de", "tables/de_results.tsv").format(contrast_id=FACTORIAL_EFFECT_Y),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "03_differential" / "factorial.R"),
        utils=UTILS_R
    output:
        effect=factorial("tables/effect_vs_effect_displayed.tsv"),
        profile=factorial("tables/interaction_profile_displayed.tsv"),
        expression=factorial("tables/group_expression_displayed.tsv"),
        effect_pdf=factorial("figures/effect_vs_effect.pdf"),
        effect_png=factorial("figures/effect_vs_effect.png"),
        profile_pdf=factorial("figures/interaction_profile.pdf"),
        profile_png=factorial("figures/interaction_profile.png"),
        expression_pdf=factorial("figures/group_expression.pdf"),
        expression_png=factorial("figures/group_expression.png"),
        summary=factorial("factorial_summary.json")
    log:
        factorial("logs/factorial.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--vst-expression {input.vst_expression:q} --samples {input.samples:q} "
        "--de-x {input.de_x:q} --de-y {input.de_y:q} "
        "--outdir {RESULTS}/factorial > {log:q} 2>&1"

rule family_fit:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        families=FAMILIES,
        estimands=ESTIMANDS,
        cell_weights=CELL_WEIGHTS,
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "03_differential" / "family_fit.R"),
        utils=UTILS_R,
        estimands_r=ESTIMANDS_R
    output:
        dds=family("objects/deseq2.rds"),
        covariance=family("tables/coefficient_covariance.tsv"),
        contrast_matrix=family("tables/contrast_matrix.tsv"),
        universe=family("tables/tested_gene_universe.tsv"),
        diagnostics=family("tables/design_diagnostics.tsv"),
        diagnostics_pdf=family("figures/design_diagnostics.pdf"),
        diagnostics_png=family("figures/design_diagnostics.png"),
        summary=family("family_summary.json")
    log:
        family("logs/family_fit.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--counts {input.counts:q} --samples {input.samples:q} "
        "--families {input.families:q} --estimands {input.estimands:q} "
        "--cell-weights {input.cell_weights:q} "
        "--family-id {wildcards.family_id:q} "
        "--outdir {RESULTS}/families/{wildcards.family_id} > {log:q} 2>&1"

rule family_estimand:
    input:
        counts=COUNTS,
        samples=SAMPLES,
        annotation=ANNOTATION,
        families=FAMILIES,
        estimands=ESTIMANDS,
        cell_weights=CELL_WEIGHTS,
        dds=lambda wildcards: str(
            RESULTS / "families" / FAMILY_OF_ESTIMAND[wildcards.contrast_id] / "objects" / "deseq2.rds"
        ),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "03_differential" / "estimand.R"),
        render=str(WORKFLOW_ROOT / "scripts" / "03_differential" / "de_render.R"),
        utils=UTILS_R,
        estimands_r=ESTIMANDS_R
    output:
        DE_PATTERNS
    log:
        analysis("de", "logs/estimand.log")
    params:
        family_dir=lambda wildcards: str(
            RESULTS / "families" / FAMILY_OF_ESTIMAND[wildcards.contrast_id]
        )
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--counts {input.counts:q} --samples {input.samples:q} "
        "--annotation {input.annotation:q} --families {input.families:q} "
        "--estimands {input.estimands:q} --cell-weights {input.cell_weights:q} "
        "--family-dir {params.family_dir:q} --estimand-id {wildcards.contrast_id:q} "
        "--outdir {RESULTS}/contrasts/{wildcards.contrast_id}/analyses/de > {log:q} 2>&1"

rule family_term_test:
    input:
        annotation=ANNOTATION,
        families=FAMILIES,
        term_tests=TERM_TESTS,
        dds=family("objects/deseq2.rds"),
        universe=family("tables/tested_gene_universe.tsv"),
        config=str(CONFIG_PATH),
        script=str(WORKFLOW_ROOT / "scripts" / "03_differential" / "term_test.R"),
        utils=UTILS_R,
        estimands_r=ESTIMANDS_R
    output:
        results=term_test("tables/lrt_results.tsv"),
        pvalue_table=term_test("tables/pvalue_distribution_displayed.tsv"),
        pvalue_pdf=term_test("figures/pvalue_distribution.pdf"),
        pvalue_png=term_test("figures/pvalue_distribution.png"),
        summary=term_test("term_test_summary.json")
    log:
        term_test("logs/term_test.log")
    conda:
        R_ENV
    shell:
        "Rscript --vanilla {input.script} --project-config {input.config:q} "
        "--annotation {input.annotation:q} --families {input.families:q} "
        "--term-tests {input.term_tests:q} "
        "--family-dir {RESULTS}/families/{wildcards.family_id} "
        "--family-id {wildcards.family_id:q} --term-test-id {wildcards.term_test_id:q} "
        "--outdir {RESULTS}/families/{wildcards.family_id}/term_tests/{wildcards.term_test_id} "
        "> {log:q} 2>&1"
