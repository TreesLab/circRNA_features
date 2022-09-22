
# ----- check if donor/acceptor site is at an annotated exon boundary ----- #

rule circmimi_annotation_db:
    input:
        config['annotation']
    output:
        "results/indexes/circmimi_refs/human/annotation.db"
    log:
        "results/logs/circmimi_annotation_db/log"
    conda:
        "../envs/circmimi.yaml"
    shell:
        "circmimi_tools gendb {input} {output}"


rule circRNAs_without_titles:
    input:
        config['circRNAs']
    output:
        "results/circRNAs_related/circRNAs_without_titles/circRNAs.no_title.tsv"
    log:
        "results/logs/circRNAs_without_titles/log"
    shell:
        "(cat {input} | sed '1d' > {output}) 2> {log}"


rule check_donor_acceptor_at_annotated_exon_boundary:
    input:
        circRNAs = "results/circRNAs_related/circRNAs_without_titles/circRNAs.no_title.tsv",
        anno_db = "results/indexes/circmimi_refs/human/annotation.db"
    output:
        "results/circRNAs_related/check_donor_acceptor_at_annotated_exon_boundary/circRNAs.check_annotation.tsv.raw",
        "results/circRNAs_related/check_donor_acceptor_at_annotated_exon_boundary/circRNAs.check_annotation.tsv"
    log:
        "results/logs/check_donor_acceptor_at_annotated_exon_boundary/log"
    conda:
        "../envs/circmimi.yaml"
    shell:
        "(circmimi_tools check annotation {input.anno_db} {input.circRNAs} {output[0]} "
        "  && cat {output[0]} "
        "      | awk 'BEGIN{{FS=\"\t\";OFS=\"\t\"}}{{if(($7==\"\")&&($8==\"\")&&($9==\"\")){{print $1, $2, $3, $4, $5, $6, 0, 0, 0}} else {{print $0}}}}'"
        "      > {output[1]})" # wait for hotfix of circmimi
        " 2> {log}"


# ----- check if there are alternative splicing events on the donor/acceptor site ----- #

rule check_AS_events_of_sites:
    input:
        gtf = "results/indexes/bgzipped_GTF/annotaion.gtf.bgz",
        sites = "results/all_circRNAs_sites/circRNAs.sites.tsv"
    output:
        "results/circRNAs_related/check_AS_events_of_sites/circRNAs.sites.check_AS.tsv"
    log:
        "results/logs/check_AS_events_of_sites/log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "python workflow/utils/check_AS_event.py -g {input.gtf} {input.sites} --log_file /dev/stderr > {output} 2> {log}"


rule create_table_with_AS_check:
    input:
        config['circRNAs'],
        "results/circRNAs_related/check_AS_events_of_sites/circRNAs.sites.check_AS.tsv"
    output:
        "results/circRNAs_related/table_with_AS_check/circRNAs.AS_check.tsv"
    log:
        "results/logs/create_table_with_AS_check/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/create_table_with_AS_check.py"


# ----- calculate the mean conservation scores of the nearby regions of the donor/acceptor sites ----- #

rule get_nearby_regions:
    input:
        config['circRNAs']
    output:
        "results/circRNAs_related/get_nearby_regions/circRNAs.nearby_regions.tsv"
    log:
        "results/logs/get_nearby_regions/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/get_nearby_regions.py"


rule get_mean_consevation_scores:
    input:
        "results/circRNAs_related/get_nearby_regions/circRNAs.nearby_regions.tsv",
        config['phyloP'],
        config['phastCons']
    output:
        "results/circRNAs_related/get_mean_consevation_scores/circRNAs.nearby_regions.mean_phyloP_phastCons.tsv"
    log:
        "results/logs/get_mean_consevation_scores/log"
    conda:
        "../envs/python3.yaml"
    shell:
        "(cat {input[0]} "
        "  | python workflow/utils/conservation_analysis/get_mean_phyloP.py -bw {input[1]} - "
        "  | python workflow/utils/conservation_analysis/get_mean_phyloP.py -bw {input[2]} - "
        "  > {output}) "
        " 2> {log}"


rule create_table_with_conservation_scores:
    input:
        "results/circRNAs_related/get_mean_consevation_scores/circRNAs.nearby_regions.mean_phyloP_phastCons.tsv"
    output:
        "results/circRNAs_related/table_with_conservation_scores/circRNAs.conservation_scores.tsv"
    log:
        "results/logs/create_table_with_conservation_scores/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/create_table_with_conservation_scores.py"


# ----- calculate the splicing scores of the donor/acceptor sites ----- #

rule get_bed_for_MaxEntScan:
    input:
        config['circRNAs']
    output:
        five_prime_bed = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.5p.MaxEntScan.bed",
        three_prime_bed = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.3p.MaxEntScan.bed"
    log:
        "results/logs/get_bed_for_MaxEntScan/log"
    conda:
        "../envs/python3.yaml"
    script:
        "../scripts/get_bed_for_MaxEntScan.py"

rule get_fasta_for_MaxEntScan:
    input:
        genome = config['genome'],
        five_prime_bed = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.5p.MaxEntScan.bed",
        three_prime_bed = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.3p.MaxEntScan.bed"
    output:
        five_prime_fa = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.5p.MaxEntScan.fa",
        three_prime_fa = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.3p.MaxEntScan.fa"
    log:
        "results/logs/get_fasta_for_MaxEntScan/log"
    conda:
        "../envs/mapping.yaml"
    shell:
        "(bedtools getfasta -fi {input.genome} -bed {input.five_prime_bed} -s -name -fo {output.five_prime_fa} &&"
        "  bedtools getfasta -fi {input.genome} -bed {input.three_prime_bed} -s -name -fo {output.three_prime_fa})"
        " 2> {log}"

# rule predict_splicing_scores_by_MaxEntScan:
#     input:
#         five_prime_fa = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.5p.MaxEntScan.fa",
#         three_prime_fa = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.3p.MaxEntScan.fa"
#     output:
#         five_prime_results = "results/circRNAs_related/predict_splicing_scores_by_MaxEntScan/circRNAs.5p.MaxEntScan_results",
#         three_prime_results = "results/circRNAs_related/predict_splicing_scores_by_MaxEntScan/circRNAs.3p.MaxEntScan_results"
#     log:
#         "results/logs/predict_splicing_scores_by_MaxEntScan/log"
#     conda:
#         "../envs/MaxEntScan.yaml"
#     shell:
#         "(maxentscan_score5.pl {input.five_prime_fa} > {output.five_prime_results} && "
#         "  maxentscan_score3.pl {input.three_prime_fa} > {output.three_prime_results})"
#         " 2> {log}"


rule predict_splicing_scores_by_MaxEntScan_web_server:
    input:
        five_prime_fa = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.5p.MaxEntScan.fa",
        three_prime_fa = "results/circRNAs_related/fasta_for_MaxEntScan/circRNAs.3p.MaxEntScan.fa"
    output:
        five_prime_results = "results/circRNAs_related/predict_splicing_scores_by_MaxEntScan_web_server/circRNAs.5p.MaxEntScan_results",
        three_prime_results = "results/circRNAs_related/predict_splicing_scores_by_MaxEntScan_web_server/circRNAs.3p.MaxEntScan_results"
    log:
        "results/logs/predict_splicing_scores_by_MaxEntScan_web_server/log"
    conda:
        "../envs/web.yaml"
    script:
        "../scripts/predict_splicing_scores_by_MaxEntScan_web_server.py"

rule create_table_with_splicing_scores:
    input:
        five_prime_results = "results/circRNAs_related/predict_splicing_scores_by_MaxEntScan_web_server/circRNAs.5p.MaxEntScan_results",
        three_prime_results = "results/circRNAs_related/predict_splicing_scores_by_MaxEntScan_web_server/circRNAs.3p.MaxEntScan_results"
    output:
        "results/circRNAs_related/table_with_splicing_scores/circRNAs.MaxEntScan_scores.tsv"
    log:
        "results/logs/create_table_with_splicing_scores/log"
    conda:
        "../envs/python3.yaml"
    shell:
        "(echo \"event_id\tMAXENT(donor)\tMM(donor)\tWMM(donor)\tMAXENT(acceptor)\tMM(acceptor)\tWMM(acceptor)\";"
        " paste {input.five_prime_results} {input.three_prime_results} | cut -f '1,3-5,8-') > {output} 2> {log}"

