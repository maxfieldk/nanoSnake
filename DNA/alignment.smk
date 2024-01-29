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
        targetAlignment = "intermediates/{sample}/alignments/{sample}_to_{target}.bam",
        targetAlignmentindex = "intermediates/{sample}/alignments/{sample}_to_{target}.bam.bai",
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
cat {params.refGenome} {params.targetFasta} > temp_{wildcards.sample}_{wildcards.target}.fa
minimap2 -ax map-ont -t 20 temp_{wildcards.sample}_{wildcards.target}.fa {input.fastq} \
| samtools sort -@4 -T $(dirname {output.targetAlignment}) -O bam -o {output.targetAlignment}
samtools index -@8 {output.targetAlignment}
        """

rule reAlignFqsExtractRegion:
    input:
        targetAlignment = "intermediates/{sample}/alignments/{sample}_to_{target}.bam",
    output:
        targetAlignmentROI = "intermediates/{sample}/alignments/{sample}_to_{target}.bam",
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
cat {params.refGenome} {params.targetFasta} > temp_{wildcards.sample}_{wildcards.target}.fa
minimap2 -ax map-ont -t 20 temp_{wildcards.sample}_{wildcards.target}.fa {input.fastq} \
| samtools sort -@4 -T $(dirname {output.targetAlignment}) -O bam -o {output.targetAlignment}
samtools index -@8 {output.targetAlignment}
        """

# rule reAlignFqs:
#     input:
#         fastq = "intermediates/{sample}/fastqs/{sample}.fq.gz"
#     output:
#         targetAlignment = "intermediates/{sample}/alignments/{sample}_to_{target}.bam",
#         targetAlignmentindex = "intermediates/{sample}/alignments/{sample}_to_{target}.bam.bai",
#     params:
#         targetFasta = lambda wildcards: config["{}Fasta".format(wildcards.target)]
#     conda: "omics"
#     resources:
#         cpus_per_task =32,
#         threads = 32,
#         runtime = 600,
#         mem_mb = 128000,
#     shell:
#         """
# minimap2 -ax map-ont -t 20 {params.targetFasta} {input.fastq} \
# | samtools view -b -F 2308 \
# | samtools sort -@4 -T $(dirname {output.targetAlignment}) -O bam -o {output.targetAlignment}
# samtools index -@8 {output.targetAlignment}
#         """


# rule reAlignToInsert:
#     input:
#         fastq = "intermediates/{sample}/fastqs/{sample}.fq"
#     output:
#         plasmidAlignment = "intermediates/{sample}/alignments/{sample}_to_plasmid.bam",
#         insertAlignment = "intermediates/{sample}/alignments/{sample}_to_insert.bam",
#         backboneAlignment = "intermediates/{sample}/alignments/{sample}_to_backbone.bam",
#         plasmidAlignmentindex = "intermediates/{sample}/alignments/{sample}_to_plasmid.bam.bai",
#         insertAlignmentindex = "intermediates/{sample}/alignments/{sample}_to_insert.bam.bai",
#         backboneAlignmentindex = "intermediates/{sample}/alignments/{sample}_to_backbone.bam.bai"
#     params:
#         plasmidFasta = config["plasmidFasta"],
#         insertFasta = config["insertFasta"],
#         backboneFasta = config["backboneFasta"]
#     conda: "omics"
#     resources:
#         cpus_per_task =32,
#         threads = 32,
#         runtime = 600,
#         mem_mb = 128000,
#     shell:
#         """
# minimap2 -ax map-ont -t 20 {params.plasmidFasta} {input.fastq} \
# | samtools sort -@4 -T $(dirname {output.plasmidAlignment}) -O bam -o {output.plasmidAlignment}
# samtools index -@8 {output.plasmidAlignment}

# minimap2 -ax map-ont -t 20 {params.insertFasta} {input.fastq} \
# | samtools sort -@4 -T $(dirname {output.insertAlignment}) -O bam -o {output.insertAlignment}
# samtools index -@8 {output.insertAlignment}

# minimap2 -ax map-ont -t 20 {params.backboneFasta} {input.fastq} \
# | samtools sort -@4 -T $(dirname {output.backboneAlignment}) -O bam -o {output.backboneAlignment}
# samtools index -@8 {output.backboneAlignment}
#         """

rule extractReadIDs:
    input:
        bam = "{path}.bam",
    output:
        readIDs = "{path}.readIDs.txt",
    conda: "omics"
    resources:
        cpus_per_task =10,
        threads = 10,
        runtime = 600,
        mem_mb = 64000,
    shell:
        """
samtools view {input.bam} | awk -F'\t' '{{print $1}}' > {output.readIDs}
        """

# rule extractMappedAlignments:
#     input:
#         bam = "{path}.bam",
#     output:
#         mappedbam = "{path}.mapped.bam",
#     conda: "omics"
#     resources:
#         cpus_per_task =10,
#         threads = 10,
#         runtime = 600,
#         mem_mb = 64000,
#     shell:
#         """
# samtools view -b -F 4 -F 0x100 -F 0x800 {input.bam} > {output.mappedbam}
#         """

rule subsetBamByReadIDs:
    input:
        genomebam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        readIDs = "intermediates/{sample}/alignments/{sample}_to_{target}.readIDs.txt"
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





