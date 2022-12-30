
rule simplify_circRNAs_format:
    input:
        config['circRNAs']
    output:
        "results/circRNAs_related/simplify_circRNAs_format/circRNAs.tsv"
    log:
        "results/logs/simplify_circRNAs_format/log"
    shell:
        "(cat {input} | sed '1d' | cut -f '-4' > {output}) 2> {log}"


rule get_RCS:
    input:
        config['genome'],
        "results/circRNAs_related/simplify_circRNAs_format/circRNAs.tsv"
    output:
        "results/circRNAs_related/get_RCS/RCS.tsv"
    log:
        "results/logs/get_RCS/log"
    params:
        dist = 20000
    threads: 30
    conda:
        "../envs/mapping.yaml"
    shell:
        "python workflow/utils/get_RCS/get_RCS.py {input} --dist {params.dist} -p {threads} > {output} 2> {log}"


rule RCS_filter:
    input:
        "results/circRNAs_related/get_RCS/RCS.tsv"
    output:
        "results/circRNAs_related/RCS_filter/RCS.filtered.tsv"
    log:
        "results/logs/RCS_filter/log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "python workflow/utils/get_RCS/RCS_filter.py {input} > {output} 2> {log}"


rule get_RCS_summary:
    input:
        "results/circRNAs_related/RCS_filter/RCS.filtered.tsv"
    output:
        "results/circRNAs_related/get_RCS_summary/RCS.filtered.summary.tsv"
    log:
        "results/logs/get_RCS_summary/log"
    conda:
        "../envs/mapping.yaml"
    params:
        dist = 20000
    shell:
        "(python workflow/utils/get_RCS/get_RCS_summary.py {input} --dist {params.dist} > {output}) 2> {log}"


rule get_RCS_status:
    input:
        "results/circRNAs_related/get_RCS_summary/RCS.filtered.summary.tsv"




rule get_flanking_bed:
    input:
        config['circRNAs']
    output:
        "results/circRNAs_related/flanking_bed/circRNAs.flanking_region.bed"
    log:
        "results/logs/get_flanking_bed/log"
    conda:
        "../envs/python3.yaml"
    params:
        dist = 1000
    script:
        "../scripts/get_flanking_bed.py"


rule get_RBP_on_flanking_1kb:
    input:
        "results/circRNAs_related/flanking_bed/circRNAs.flanking_region.bed",
        config['ENCORI_RBP']
    output:
        "results/circRNAs_related/get_RBP_on_flanking_1kb/circRNAs.flanking_region.bed.intersect_RBP"
    log:
        "results/logs/get_RBP_on_flanking_1kb/log"
    conda:
        "../envs/circmimi.yaml"
    shell:
        "(bedtools intersect -wao -a {input[0]} -b {input[1]} > {output}) 2> {log}"
    

rule calc_RBP_coverage:
    input:
        "results/circRNAs_related/get_RBP_on_flanking_1kb/circRNAs.flanking_region.bed.intersect_RBP"
    output:
        "results/circRNAs_related/calc_RBP_coverage/circRNAs.flanking_region.bed.intersect_RBP.coverage"
    log:
        "results/logs/calc_RBP_coverage/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/calc_RBP_coverage.py"


rule get_common_RBP:
    input:
        "results/circRNAs_related/calc_RBP_coverage/circRNAs.flanking_region.bed.intersect_RBP.coverage"
    output:
        "results/circRNAs_related/get_common_RBP/circRNAs.flanking_region.bed.intersect_RBP.coverage.common_RBP"
    log:
        "results/logs/get_common_RBP/log"
    conda:
        "../envs/python3.yaml"
    shell:
        "(cat {input} "
        "  | python workflow/utils/flanking_region_RBP/simplify_format.py "
        "  | python workflow/utils/flanking_region_RBP/get_common_RBP.py "
        "  > {output}) "
        " 2> {log}"


rule RBP_pairs_on_flanking_1kb:
    input:
        "results/circRNAs_related/get_common_RBP/circRNAs.flanking_region.bed.intersect_RBP.coverage.common_RBP"


rule get_RBP_dist_data:
    input:
        "results/circRNAs_related/calc_RBP_coverage/circRNAs.flanking_region.bed.intersect_RBP.coverage"
    output:
        "results/circRNAs_related/get_RBP_dist_data/circRNAs.flanking_region.bed.intersect_RBP.coverage.dist_data"
    log:
        "results/logs/get_RBP_dist_data/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/get_RBP_dist_data.py"

rule get_common_RBP_with_dist_data:
    input:
        "results/circRNAs_related/get_RBP_dist_data/circRNAs.flanking_region.bed.intersect_RBP.coverage.dist_data"
    output:
        "results/circRNAs_related/get_RBP_dist_data/circRNAs.flanking_region.bed.intersect_RBP.coverage.dist_data.common_RBP"
    log:
        "results/logs/get_common_RBP_with_dist_data/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/get_common_RBP_with_dist_data.py"

rule get_min_dist_RBP:
    input:
        "results/circRNAs_related/get_RBP_dist_data/circRNAs.flanking_region.bed.intersect_RBP.coverage.dist_data.common_RBP"
    output:
        "results/circRNAs_related/get_RBP_dist_data/circRNAs.flanking_region.bed.intersect_RBP.coverage.dist_data.common_RBP.min_dist_RBP"
    log:
        "results/logs/get_min_dist_RBP/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/get_min_dist_RBP.py"
