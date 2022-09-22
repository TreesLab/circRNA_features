

# ----- append the table from the circAtlas web site ----- #

rule create_table_from_circAtlas:
    input:
        config['circAtlas']
    output:
        "results/circRNAs_related/circAtlas_table/circAtlas_table.human.2.tsv"
    log:
        "results/logs/create_table_from_circAtlas/log"
    conda:
        "../envs/python3.yaml"
    shell:
        "(python workflow/utils/generate_circAtlas_table.py {input} > {output}) 2> {log}"


# ----- check the database transCirc for the circRNAs-junction supporting evidences ----- #

rule create_table_from_transCirc:
    input:
        config['transCirc'],
        "results/circRNAs_related/circAtlas_table/circAtlas_table.human.2.tsv"
    output:
        "results/circRNAs_related/transCirc_table/transCirc_table.tsv"
    log:
        "results/logs/create_table_from_transCirc/log"
    conda:
        "../envs/python3.yaml"
    script:
        "scripts/create_table_from_transCirc.py"


