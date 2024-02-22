rule extractFastqs:
    input:
        bam = "intermediates/{sample}/alignments/{sample}.sorted.bam"
    output:
        fastq = "intermediates/{sample}/fastqs/{sample}.fq"
    resources:
        cpus_per_task =10,
        threads = 10,
        runtime = 300,
        mem_mb = 32000
    conda: "omics"
    shell:
        """
samtools fastq -c6 -@8 {input.bam} > {output.fastq}
        """

rule reAlignFqs:
    input:
        fastq = "intermediates/{sample}/fastqs/{sample}.fq"
    output:
        targetAlignment = "intermediates/{sample}/alignments/{sample}_to_{target, [A-Za-z0-9]+}.bam",
        targetAlignmentindex = "intermediates/{sample}/alignments/{sample}_to_{target, [A-Za-z0-9]+}.bam.bai",
    params:
        refGenome = config["reference"],
        targetFasta = lambda wildcards: config["{}Fasta".format(wildcards.target)]
    conda: "omics"
    resources:
        cpus_per_task =32,
        threads = 32,
        runtime = 1200,
        mem_mb = 128000,
    shell:
        """
cat {params.refGenome} {params.targetFasta} > temp_{wildcards.sample}_{wildcards.target}.fa
minimap2 -ax map-ont -t 20 temp_{wildcards.sample}_{wildcards.target}.fa {input.fastq} \
| samtools sort -@4 -T $(dirname {output.targetAlignment}) -O bam -o {output.targetAlignment}
samtools index -@8 {output.targetAlignment}
        """

rule reAlignFqsExtractRegion:
    input:
        targetAlignment = "intermediates/{sample}/alignments/{sample}_to_{target}.bam",
    output:
        targetAlignmentROI = "intermediates/{sample}/alignments/{sample}_to_{target}_ROI.bam",
    params:
        refGenome = config["reference"],
        targetFasta = lambda wildcards: config["{}Fasta".format(wildcards.target)]
    conda: "omics"
    resources:
        cpus_per_task =32,
        threads = 32,
        runtime = 600,
        mem_mb = 128000,
    shell:
        """
samtools view -h -b {input.targetAlignment} "$(basename {params.targetFasta})" > {output.targetAlignmentROI}
samtools index {output.targetAlignmentROI}
        """


rule extractReadIDs:
    input:
        targetAlignmentROI = "intermediates/{sample}/alignments/{sample}_to_{target}_ROI.bam"
    output:
        readIDs = "intermediates/{sample}/alignments/{sample}_to_{target}_ROI.readIDs.txt"
    conda: "omics"
    resources:
        cpus_per_task =10,
        threads = 10,
        runtime = 600,
        mem_mb = 64000,
    shell:
        """
samtools view {input.targetAlignmentROI} | awk -F'\t' '{{print $1}}' > {output.readIDs}
        """


rule subsetBamByReadIDs:
    input:
        genomebam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        readIDs = "intermediates/{sample}/alignments/{sample}_to_{target}_ROI.readIDs.txt"
    output:
        bam = "intermediates/{sample}/alignments/{sample}_genome_extracted_reads_that_mapped_to_{target}.bam",
    conda: "omics"
    resources:
        cpus_per_task =10,
        threads = 10,
        runtime = 600,
        mem_mb = 64000,
    shell:
        """
samtools view -b -F 2304 -N {input.readIDs} {input.genomebam} > {output.bam}
        """





