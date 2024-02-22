rule sortBam:
    input:
        bam = "{bampath}.bam"
    output:
        sortedbam =  "{bampath}.sorted.bam",
    resources:
        cpus_per_task =10,
        mem_mb = 64000
    conda:
        "omics"
    shell:
        """
samtools sort -@8 -m4g {input.bam} > {output.sortedbam}
        """

rule IndexBam:
    input:
        bam = "{bampath}.bam"
    output:
        index = "{bampath}.bam.bai"
    resources:
        cpus_per_task =10,
        mem_mb = 64000
    conda:
        "omics"
    shell:
        """
samtools index  -@6 {input.bam}
        """

# rule common_samtobam:
#     input:
#         sam = "{path}.sam"
#     output:
#         bam =  "{path}.bam",
#     resources:
#         cpus_per_task =10,
#         mem_mb = 30000
#     conda:
#         "omics"
#     shell:
#         """
# samtools view -@8 {input.sam} > {output.bam}
#         """

#filter out secondary, keep primary and supplementary, and do some quality score filtering
rule filterbam:
    input:
        bam = "{path}.sorted.bam"
    output:
        bam = "{path}.sorted.filtered.bam"
    resources:
        cpus_per_task =6,
        mem_mb = 40000
    conda:
        "omics"
    shell: "samtools view -b -F 0x100 -q 10 -e '[qs] > 10' {input.bam} > {output.bam}"

rule mergeBams:
    input:
        expand("intermediates/{sample}/alignments/{sample}.sorted.bam", sample = samples)
    output:
        mergedbam = "intermediates/merged.sorted.bam"
    conda:
        "minimap2"
    resources:
        cpus_per_task =8,
        mem_mb = 128000,
        slurm_extra="--time=2:00:00 --constraint=cascade"
    log: "logs/mergeBams.log"
    shell: "samtools merge --threads 6 {output.mergedbam} {input} 2> {log}"
    

