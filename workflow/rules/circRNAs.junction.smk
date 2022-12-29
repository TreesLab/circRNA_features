


# ----- check if there are ambiguous alignments for the circRNA sequences nearby junction ----- #

rule generate_other_transcripts_reference:
    input:
        config['genome'],
        config['pc_ref'],
        config['lncRNA_ref']
    output:
        "results/indexes/circmimi_refs/human/other.fa"
    log:
        "results/logs/generate_other_transcripts_reference/log"
    conda:
        "../envs/python3.yaml"
    params:
        src_type = 'ensembl'
    shell:
        "python workflow/utils/amb_check/generate_other_ref.py {input} {output} --src_type {params.src_type} 2> {log}"
        

rule split_pseudo_reference:
    input:
        "results/pseudo_reference/pseudo_reference.{L}.fa"
    output:
        directory("results/circRNAs_related/split_pseudo_reference/pseudo_reference.{L}.fa.parts.each_{num_records}/")
    log:
        "results/logs/split_pseudo_reference/{L}.{num_records}.log"
    params:
        lines = lambda wildcards: 2 * int(wildcards.num_records),
        suffix_length = 5,
        out_prefix = lambda wildcards, output: f"{output}/pseudo_reference.{wildcards.L}.part_"
    shell:
        "(mkdir -p {output} && "
        "  cat {input} | split -l {params.lines} --numeric-suffixes=1 --additional-suffix='.fa' -a {params.suffix_length} - {params.out_prefix})"
        " 2> {log}"


rule check_ambiguous_alignments:
    input:
        config['genome'],
        "results/indexes/circmimi_refs/human/other.fa",
        "results/circRNAs_related/split_pseudo_reference/pseudo_reference.100.fa.parts.each_100/"
    output:
        directory("results/circRNAs_related/check_ambiguous_alignments/pseudo_reference.100.fa.parts.results/")
    log:
        "results/logs/check_ambiguous_alignments/log"
    conda:
        "../envs/circmimi.yaml"
    threads: 30
    params:
        checkAA = "workflow/utils/amb_check/checkAA_reads.py",
        run_all_checkAA = "workflow/utils/amb_check/run_all_checkAA.py"
    shell:
        "(ls {input[2]}/* "
        "  | python {params.run_all_checkAA} -rg {input[0]} -ro {input[1]} - {output} -p {threads} --checkAA_bin {params.checkAA}) "
        " 2> {log}"


rule merge_checkAA_results:
    input:
        "results/circRNAs_related/check_ambiguous_alignments/pseudo_reference.100.fa.parts.results/"
    output:
        "results/circRNAs_related/check_ambiguous_alignments/circRNAs.checkAA.tsv"
    log:
        "results/logs/merge_checkAA_results/log"
    shell:
        "cat {input}/* | head -1 > {output} && ls {input}/* | while read file; do; cat $file | sed '1d' >> {output}"



# ----- check if there are miRNA-binding sites across the circRNA junction ----- #

rule predict_miRNA_binding_sites_across_junction:
    input:
        config['mir_ref'],
        "results/circRNAs_related/split_pseudo_reference/pseudo_reference.50.fa.parts.each_100/"
    output:
        protected(directory("results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.parts.results/"))
    log:
        "results/logs/predict_miRNA_binding_sites_across_junction/log"
    conda:
        "../envs/circmimi.yaml"
    params:
        filename_prefix = "pseudo_reference.50.part_"
    threads: 30
    script:
        "../scripts/predict_miRNA_binding_sites_across_junction.py"


rule parse_miranda_results:
    input:
        "results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.parts.results/"
    output:
        directory("results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.parts.results.parsed/")
    log:
        "results/logs/parse_miranda_results/log"
    conda:
        "../envs/python3.yaml"
    params:
        filename_prefix = "pseudo_reference.50.part_"
    threads: 10
    script:
        "../scripts/parse_miranda_results.py"


rule merge_all_miranda_results:
    input:
        "results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.parts.results.parsed/"
    output:
        "results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.miranda_results"
    log:
        "results/logs/merge_all_miranda_results/log"
    conda:
        "../envs/python3.yaml"
    params:
        filename_prefix = "pseudo_reference.50.part_"
    script:
        "../scripts/merge_all_miranda_results.py"


rule get_crossing_junction_miRNA:
    input:
        "results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.miranda_results"
    output:
        "results/circRNAs_related/predict_miRNA_binding_sites_across_junction/pseudo_reference.50.fa.miranda_results.crossing_junction_miRNAs"
    log:
        "results/logs/get_crossing_junction_miRNA/log"
    conda:
        "../envs/python3.yaml"
    params:
        L = 50,
        cross_junction = 5,
        score_threshold = 155
    script:
        "../scripts/get_crossing_junction_miRNA.py"


# ----- check if there are RBP-binding sites across the circRNA junction ----- #

rule predict_RBP_binding_sites_across_junction:
    input:
        "results/circRNAs_related/split_pseudo_reference/pseudo_reference.50.fa.parts.each_5000/"
    output:
        directory("results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.parts.each_5000.results")
    log:
        "results/logs/predict_RBP_binding_sites_across_junction/log"
    conda:
        "../envs/web.yaml"
    params:
        filename_prefix = "pseudo_reference.50.part_",
        job_limit = 5
    shell:
        "python workflow/utils/RBPmap/predict_RBP_binding_sites_across_junction.py "
        " {input} {output} --prefix {params.filename_prefix} --job_limit {params.job_limit}"


rule parse_RBPmap_results:
    input:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.parts.each_5000.results"
    output:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results"
    log:
        "results/logs/parse_RBPmap_results/log"
    conda:
        "../envs/python3.yaml"
    shell:
        "python workflow/utils/parse_RBPmap_results.py {input}/* > {output}"


rule get_crossing_junction_RBP:
    input:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results"
    output:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results.crossing_junction_RBPs",
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results.crossing_junction_RBPs.count"
    log:
        "results/logs/get_crossing_junction_RBP/log"
    conda:
        "../envs/python3.yaml"
    params:
        L = 50,
        cross_junction = 2
    script:
        "../scripts/get_crossing_junction_RBP.py"


## RBPmap with "high stringency" and "conservation filter"

rule predict_RBP_binding_sites_across_junction_high_cons:
    input:
        "results/circRNAs_related/split_pseudo_reference/pseudo_reference.50.fa.parts.each_5000/"
    output:
        directory("results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.parts.each_5000.results.high_cons")
    log:
        "results/logs/predict_RBP_binding_sites_across_junction_high_cons/log"
    conda:
        "../envs/web.yaml"
    params:
        filename_prefix = "pseudo_reference.50.part_",
        job_limit = 5
    shell:
        "python workflow/utils/RBPmap/predict_RBP_binding_sites_across_junction.py "
        " {input} {output} --prefix {params.filename_prefix} --job_limit {params.job_limit} "
        "--high-cons"


rule parse_RBPmap_results_high_cons:
    input:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.parts.each_5000.results.high_cons"
    output:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results.high_cons"
    log:
        "results/logs/parse_RBPmap_results_high_cons/log"
    conda:
        "../envs/python3.yaml"
    shell:
        "python workflow/utils/parse_RBPmap_results.py {input}/* > {output}"


rule get_crossing_junction_RBP_high_cons:
    input:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results.high_cons"
    output:
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results.high_cons.crossing_junction_RBPs",
        "results/circRNAs_related/predict_RBP_binding_sites_across_junction/pseudo_reference.50.fa.RBPmap_results.high_cons.crossing_junction_RBPs.count"
    log:
        "results/logs/get_crossing_junction_RBP_high_cons/log"
    conda:
        "../envs/python3.yaml"
    params:
        L = 50,
        cross_junction = 2
    script:
        "../scripts/get_crossing_junction_RBP.py"


# ----- check if there are G-quadruplex structure across the circRNA junction ----- #

rule predict_G_quadruplex_structure_across_junction:
    input:
        "results/pseudo_reference/pseudo_reference.50.fa"
    output:
        "results/circRNAs_related/predict_G_quadruplex_structure_across_junction/QGRS_results.txt"
    log:
        "results/logs/predict_G_quadruplex_structure_across_junction/log"
    conda:
        "../envs/qgrs.yaml"
    script:
        "../scripts/predict_G_quadruplex_structure_across_junction.py"


rule get_cross_junction_G_quadruplex:
    input:
        "results/circRNAs_related/predict_G_quadruplex_structure_across_junction/QGRS_results.txt"
    output:
        "results/circRNAs_related/predict_G_quadruplex_structure_across_junction/QGRS_results.cross_junction_{dist}.txt"
    log:
        "results/logs/get_cross_junction_G_quadruplex/{dist}.log"
    conda:
        "../envs/python3.yaml"
    params:
        L = 50
    script:
        "../scripts/get_cross_junction_G_quadruplex.py"


