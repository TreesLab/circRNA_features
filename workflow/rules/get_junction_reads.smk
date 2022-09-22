import re
import os.path



rule unzip_gz_file:
    input:
        "resources/datasets/{dataset_id}/{sample_id}/{filename}.gz"
    output:
        # temp("results/unzip_gz_file/{dataset_id}/{sample_id}/{filename}")
        "results/unzip_gz_file/{dataset_id}/{sample_id}/{filename}"
    log:
        "results/logs/unzip_gz_file/{dataset_id}/{sample_id}/{filename}.log"
    resources:
        IO_limit = 1
    shell:
        "(zcat {input} > {output}) 2> {log}"


def get_unzipped_read_files(wildcard):
    sample = list(
        filter(
            lambda sample: (sample['dataset_id'] == wildcard.dataset_id) and \
                (sample['sample_id'] == wildcard.sample_id),
            SAMPLES_LIST
        )
    )[0]

    filenames = re.split(r'[|,]', sample['filenames'])
    filenames = [
        os.path.join(
            f'results/unzip_gz_file/{wildcard.dataset_id}/{wildcard.sample_id}/',
            filename.rstrip('.gz')
        )
        for filename in filenames
    ]

    return filenames


rule fastuniq_file_list:
    input:
        get_unzipped_read_files
    output:
        temp("results/fastuniq_file_list/{dataset_id}/{sample_id}/file_list.txt")
        # "results/fastuniq_file_list/{dataset_id}/{sample_id}/file_list.txt"
    log:
        "results/logs/fastuniq_file_list/{dataset_id}/{sample_id}.log"
    shell:
        "(echo {input} | tr ' ' '\n' > {output}) 2> {log}"


rule fastuniq_reads:
    input:
        "results/fastuniq_file_list/{dataset_id}/{sample_id}/file_list.txt",
        get_unzipped_read_files # hook
    output:
        # temp("results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_1.uniq.fq"),
        # temp("results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_2.uniq.fq")
        "results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_1.uniq.fq",
        "results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_2.uniq.fq"
    log:
        "results/logs/fastuniq_reads/{dataset_id}/{sample_id}.log"
    conda:
        "../envs/mapping.yaml"
    resources:
        mem_mb = 400000,
        IO_limit = 1
    shell:
        "(fastuniq -i {input[0]} -o {output[0]} -p {output[1]}) 2> {log}"


rule compress_fastuniq_results:
    input:
        "results/fastuniq_reads/{dataset_id}/{sample_id}/{filename}.uniq.fq"
    output:
        protected("results/fastuniq_reads/{dataset_id}/{sample_id}/{filename}.uniq.fq.gz")
        # "results/fastuniq_reads/{dataset_id}/{sample_id}/{filename}.uniq.fq.gz"
    log:
        "results/logs/compress_fastuniq_results/{dataset_id}/{sample_id}/{filename}.log"
    resources:
        IO_limit = 1
    shell:
        "(gzip -c {input} > {output}) 2> {log}"


rule run_all_fastuniq:
    input:
        expand(
            outputs_with_sample_info(
                "results/fastuniq_reads/{}/{}/{}_{{mate}}.uniq.fq.gz",
                SAMPLES_LIST,
                ['dataset_id', 'sample_id', 'sample_id']
            ),
            mate=[1, 2]
        )


# rule fastuniq_reads:
#     input:
#         "resources/datasets/{dataset_id}/{sample}/{sample}_f1.fq.gz",
#         "resources/datasets/{dataset_id}/{sample}/{sample}_r2.fq.gz"
#     output:
#         "results/fastuniq_reads/{dataset_id}/{sample}/{sample}_f1.uniq.fq.gz",
#         "results/fastuniq_reads/{dataset_id}/{sample}/{sample}_r2.uniq.fq.gz"
#     conda:
#         "../envs/mapping.yaml"
#     log:
#         "results/logs/fastuniq_reads/{dataset_id}/{sample}.log"
#     resources:
#         mem_mb = 200000
#     shell:
#         "python workflow/utils/run_fastuniq.py {input} {output} 2> {log}"


rule SE_merge_replicates_reads:
    input:
        get_unzipped_read_files
    output:
        "results/SE_merge_replicates_reads/{dataset_id}/{sample_id}/{sample_id}.fastq.gz"
    log:
        "results/logs/SE_merge_replicates_reads/{dataset_id}/{sample_id}/log"
    shell:
        "(cat {input} | gzip > {output}) 2> {log}"


rule STAR_generate_index:
    input:
        genome = config['genome'],
        annotation = config['annotation']
    output:
        directory("results/indexes/STAR_index_human/")
    conda:
        "../envs/mapping.yaml"
    threads: 30
    log:
        "results/logs/STAR_generate_index/STAR_index_human.log"
    shell:
        "STAR "
        "  --runThreadN {threads}"
        "  --runMode genomeGenerate"
        "  --genomeDir {output}"
        "  --genomeFastaFiles {input.genome}"
        "  --sjdbGTFfile {input.annotation}"
        " > {log}"


def get_STAR_mapping_reads(wildcard):
    sample = list(
        filter(
            lambda sample: (sample['dataset_id'] == wildcard.dataset_id) and \
                (sample['sample_id'] == wildcard.sample_id),
            SAMPLES_LIST
        )
    )[0]

    if sample['layout'] == 'single':
        fastqs = [
            "results/SE_merge_replicates_reads/{dataset_id}/{sample_id}/{sample_id}.fastq.gz"
        ]

    elif sample['layout'] == 'paired':
        fastqs = [
            "results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_1.uniq.fq.gz",
            "results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_2.uniq.fq.gz"
        ]

    return fastqs


rule STAR_mapping:
    input:
        genomeDir = "results/indexes/STAR_index_human/",
        # fastqs = [
        #     "results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_1.uniq.fq.gz",
        #     "results/fastuniq_reads/{dataset_id}/{sample_id}/{sample_id}_2.uniq.fq.gz"
        # ]
        fastqs = get_STAR_mapping_reads
    output:
        # protected(directory("results/STAR_mapping/{dataset_id}/{sample_id}/")),
        directory("results/STAR_mapping/{dataset_id}/{sample_id}/"),
        "results/STAR_mapping/{dataset_id}/{sample_id}/Aligned.out.bam",
        "results/STAR_mapping/{dataset_id}/{sample_id}/SJ.out.tab",
        "results/STAR_mapping/{dataset_id}/{sample_id}/Unmapped.out.mate1",
        "results/STAR_mapping/{dataset_id}/{sample_id}/Unmapped.out.mate2"
    threads: 30
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/STAR_mapping/{dataset_id}/{sample_id}.log"
    shell:
        "STAR "
        "  --runThreadN {threads}"
        "  --genomeDir {input.genomeDir}"
        "  --readFilesIn {input.fastqs}"
        "  --readFilesCommand zcat"
        "  --outFileNamePrefix {output[0]}/"
        "  --outSAMtype BAM Unsorted"
        "  --outSAMunmapped Within KeepPairs"
        "  --outReadsUnmapped Fastx"
        "  --quantMode TranscriptomeSAM GeneCounts"
        "  --outSAMattributes All"
        " > {log} "
        "&& touch {output[2]}"


rule run_all_STAR_mapping:
    input:
        outputs_with_sample_info(
            "results/STAR_mapping/{}/{}/",
            SAMPLES_LIST,
            ['dataset_id', 'sample_id']
        )


rule sort_STAR_unmapped_by_read_name:
    input:
        "results/STAR_mapping/{dataset_id}/{sample_id}/Unmapped.out.mate{mate}"
    output:
        "results/STAR_mapping/{dataset_id}/{sample_id}/Unmapped.out.sorted.mate{mate}"
    log:
        "results/logs/sort_STAR_unmapped_by_read_name/{dataset_id}/{sample_id}/{mate}.log"
    shell:
        "cat {input[0]} | paste - - - - | sort -k1,1 | tr '\t' '\n' > {output[0]}"



rule generate_pseudo_reference:
    input:
        config['genome'],
        config['circRNAs']
    output:
        "results/pseudo_reference/pseudo_reference.{L}.fa"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/generate_pseudo_reference/{L}.log"
    wildcard_constraints:
        L = "[1-9][0-9]*"
    shell:
        "(cat {input[1]} | sed '1d' | cut -f '-4'"
        "  | python workflow/utils/get_flanking_seq.py {input[0]} /dev/stdin {output} {wildcards.L})"
        " 2> {log}"


rule BWA_generate_index:
    input:
        "results/pseudo_reference/pseudo_reference.{L}.fa"
    output:
        "results/indexes/bwa_index_human_pseudo_reference/pseudo_reference.{L}.fa"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/BWA_generate_index/{L}.log"
    shell:
        "cp {input} {output} && bwa index {output} 2> {log}"



rule BWA_mapping:
    input:
        "results/indexes/bwa_index_human_pseudo_reference/pseudo_reference.100.fa",
        "results/STAR_mapping/{dataset_id}/{sample_id}/Unmapped.out.sorted.mate{mate}",
        "results/STAR_mapping/{dataset_id}/{sample_id}/"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.bam"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/BWA_mapping/{dataset_id}/{sample_id}.{mate}.log"
    threads: 15
    shell:
        "(bwa mem -t {threads} -O 100 {input[0]} {input[1]} | samtools view -Sbh - > {output}) 2> {log}"


rule run_all_BWA_mapping:
    input:
        expand(
            outputs_with_sample_info(
                "results/BWA_mapping/{}/{}/bwa_unmapped_single_{{mate}}/Aligned.out.bam",
                SAMPLES_LIST,
                ['dataset_id', 'sample_id']
            ),
            mate=[1, 2]
        )


rule append_Z3_XS:
    input:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.bam"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.Z3.XS.bam"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/append_Z3_XS/{dataset_id}/{sample_id}.{mate}.log"
    shell:
        "(samtools view -h {input}"
        "  | python workflow/utils/append_Z3_tag.py"
        "  | python workflow/utils/append_similarity_XS.py"
        "  | samtools view - -hb"
        "  > {output})"
        " 2> {log}"



rule get_uniq_match:
    input:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.Z3.XS.bam"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.Z3.XS.uniq_matches.bam"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/get_uniq_match/{dataset_id}/{sample_id}.{mate}.log"
    shell:
        "(samtools view -h {input}"
        "  | python workflow/utils/get_uniq_matches.py"
        "  | samtools view - -hb"
        "  > {output})"
        " 2> {log}"


rule get_junction_reads:
    input:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.Z3.XS.uniq_matches.bam"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_{mate}/Aligned.out.Z3.XS.uniq_matches.bam.junc_reads"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/get_junction_reads/{dataset_id}/{sample_id}.{mate}.log"
    shell:
        "(samtools view -h {input}"
        "  | python workflow/utils/get_junc_reads.py 100 10 20 0.8"
        "  > {output})"
        " 2> {log}"



rule merge_junction_reads:
    input:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_1/Aligned.out.Z3.XS.uniq_matches.bam.junc_reads",
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_unmapped_single_2/Aligned.out.Z3.XS.uniq_matches.bam.junc_reads"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_junc_reads/all_junc_reads_s1_s2.tsv"
    log:
        "results/logs/merge_junction_reads/{dataset_id}/{sample_id}.log"
    shell:
        "(cat"
        "   <(cat {input[0]} | awk -F'\t' '{{print $0\"\t1\"}}')"
        "   <(cat {input[1]} | awk -F'\t' '{{print $0\"\t2\"}}')"
        "  | sort -k 1,1 -k 2,2 -k 9,9"
        "  > {output})"
        " 2> {log}"


rule retain_uniq_read_ref_pairs:
    input:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_junc_reads/all_junc_reads_s1_s2.tsv"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_junc_reads/all_junc_reads_s1_s2.tsv.uniq_read_ref"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/retain_uniq_read_ref_pairs/{dataset_id}/{sample_id}.log"
    shell:
        "(cat {input} | python workflow/utils/retain_uniq_read_ref.py > {output}) 2> {log}"


rule merge_rows:
    input:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_junc_reads/all_junc_reads_s1_s2.tsv.uniq_read_ref"
    output:
        "results/BWA_mapping/{dataset_id}/{sample_id}/bwa_junc_reads/all_junc_reads.count"
    conda:
        "../envs/mapping.yaml"
    log:
        "results/logs/merge_rows/{dataset_id}/{sample_id}.log"
    shell:
        "(cat {input}"
        "  | awk 'BEGIN{{FS=\"\t\"; OFS=\"\t\"}}{{print $2,$1}}'"
        "  | sort -k1,1 -k2,2"
        "  | python workflow/utils/merge_rows.py"
        "  > {output})"
        " 2> {log}"


rule get_all_junction_reads_count:
    input:
        outputs_with_sample_info(
            "results/BWA_mapping/{}/{}/bwa_junc_reads/all_junc_reads.count",
            SAMPLES_LIST,
            ['dataset_id', 'sample_id']
        )



def junc_reads_count_tables(wildcard):
    samples = list(
        filter(
            lambda sample: (sample['dataset_id'] == wildcard.dataset_id) and \
                (sample['species'] == wildcard.species) and \
                (sample['tissue'] == wildcard.tissue),
            SAMPLES_LIST
        )
    )

    sample_order = ('totalRNA', 'RNaseR', 'polyA_minus')
    sorted_samples = sorted(samples, key=lambda sample: sample_order.index(sample['treatments']))

    tables = [
        "results/BWA_mapping/{}/{}/bwa_junc_reads/all_junc_reads.count".format(sample['dataset_id'], sample['sample_id']) 
        for sample in sorted_samples
    ]

    return tables


rule create_junction_reads_table:
    input:
        config['circRNAs'],
        junc_reads_count_tables
    output:
        "results/junction_reads_table/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.tsv"
    conda:
        "../envs/python3.yaml"
    log:
        "results/logs/create_junction_reads_table/{dataset_id}/{species}/{tissue}/log"
    params:
        column_names = get_treatment_names
    shell:
        "python workflow/utils/merge_junc_reads_result.py {input} -c {params.column_names} > {output} 2> {log}"


rule retain_circRNAs_with_totalRNA_support:
    input:
        "results/junction_reads_table/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.tsv"
    output:
        "results/junction_reads_table/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.totalRNA.tsv"
    log:
        "results/logs/retain_circRNAs_with_totalRNA_support/{dataset_id}/{species}/{tissue}/log"
    conda:
        "../envs/python3.yaml"
    params:
        threshold = 1
    script:
        "../scripts/retain_circRNAs_with_totalRNA_support.py"


rule get_all_junction_reads_table:
    input:
        [
            "results/junction_reads_table/{}/{}/{}/circRNAs.junc_reads.count.totalRNA.tsv".format(dataset_id, species, tissue)
            for dataset_id, species, tissue in ALL_SAMPLE_GROUPS
        ]


rule count_the_number_of_reads:
    input:
        get_STAR_mapping_reads
    output:
        "results/count_the_number_of_reads/{dataset_id}/{sample_id}/num_reads.txt"
    log:
        "results/logs/count_the_number_of_reads/{dataset_id}/{sample_id}/log"
    shell:
        "(zcat {input[0]} | paste - - - - | wc -l > {output}) 2> {log}"



def num_reads_of_samples(wildcard):
    samples = list(
        filter(
            lambda sample: (sample['dataset_id'] == wildcard.dataset_id) and \
                (sample['species'] == wildcard.species) and \
                (sample['tissue'] == wildcard.tissue),
            SAMPLES_LIST
        )
    )

    sample_order = ('totalRNA', 'RNaseR', 'polyA_minus')
    sorted_samples = sorted(samples, key=lambda sample: sample_order.index(sample['treatments']))

    num_reads = [
        "results/count_the_number_of_reads/{}/{}/num_reads.txt".format(sample['dataset_id'], sample['sample_id']) 
        for sample in sorted_samples
    ]

    return num_reads



rule calc_RPM:
    input:
        "results/junction_reads_table/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.totalRNA.tsv",
        num_reads_of_samples
    output:
        "results/calc_RPM/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.totalRNA.RPM.tsv"
    log:
        "results/logs/calc_RPM/{dataset_id}/{species}/{tissue}/log"
    conda:
        "../envs/python3.yaml"
    params:
        treatments = get_treatment_names
    script:
        "../scripts/calc_RPM.py"


rule calc_all_RPM:
    input:
        [
            "results/calc_RPM/{}/{}/{}/circRNAs.junc_reads.count.totalRNA.RPM.tsv".format(dataset_id, species, tissue)
            for dataset_id, species, tissue in ALL_SAMPLE_GROUPS
        ]


rule get_circRNAs_status_in_RNAseq:
    input:
        "results/calc_RPM/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.totalRNA.RPM.tsv"
    output:
        "results/circRNAs_status_in_RNAseq/{dataset_id}/{species}/{tissue}/circRNAs.junc_reads.count.totalRNA.RPM.status.tsv"
    log:
        "results/logs/get_circRNAs_status_in_RNAseq/{dataset_id}/{species}/{tissue}/log"
    conda:
        "../envs/python3.yaml"
    params:
        treatments = get_treatment_names
    script:
        "../scripts/get_circRNAs_status_in_RNAseq.py"

rule get_all_circRNAs_status_in_RNAseq:
    input:
        [
            "results/circRNAs_status_in_RNAseq/{}/{}/{}/circRNAs.junc_reads.count.totalRNA.RPM.status.tsv".format(dataset_id, species, tissue)
            for dataset_id, species, tissue in ALL_SAMPLE_GROUPS
        ]
