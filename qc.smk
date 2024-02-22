rule multiqc:
    conda:
        "omics"
    output:
        report = "qc/multiqc_report.html",
    shell:
        """
mkdir -p {params.outputdir}
multiqc --force --outdir $(dirname {output.report}) --filename {output.report} --export .
        """

rule dorado_seqsummary:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.bam.bai"
    params:
        dorado = config["dorado"]
    output:
        "qc/{sample}/{sample}.doradosummary.txt"
    conda:
        "omics"
    shell:
        """
{params.dorado} summary {input.sortedbam} > {output}
        """

rule pycoQC:
    input:
        seqsummary = "qc/{sample}/{sample}.doradosummary.txt",
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.bam",

    output:
        "qc/{sample}/{sample}pycoQC.html"
    conda:
        "pycoQC"
    shell:
        """
pycoQC --summary_file {input.seqsummary} --bam_file {input.sortedbam} --html_outfile {output} --min_pass_qual 10
        """

rule samtoolsstats:
    input:
        bam = "{bampath}.bam"
    output:
        stats = "{bampath}.bam.stats.txt"
    resources:
        cpus_per_task =10,
        mem_mb = 64000
    conda:
        "omics"
    shell:
        """
samtools stats -@8 {input.bam} > {output.stats}
        """

rule plotbamstats:
    input:
        stats = "{bampath}.bam.stats.txt"
    output:
        bamstats = touch("{bampath}.bam.stats.plots.out")
    resources:
        cpus_per_task =10,
        mem_mb = 64000
    conda:
        "omics"
    shell:
        """
plot-bamstats -p {wildcards.bampath}.bam.stats.plots {input.stats}
        """

rule flagStat:
    input:
        bam = "{path}.bam",
    output:
        stats = "{path}.bam.flagstat.txt"
    conda: "omics"
    resources:
        cpus_per_task =10,
        threads = 10,
        runtime = 600,
        mem_mb = 64000
    shell:
        """
samtools flagstat -@10 {input.bam} > {output.stats}
        """

rule seqfuStats:
    input:
        bam = "intermediates/{sample}/fastqs/{sample}.fq.gz",
    output:
        stats = "qc/{sample}_seqfuStats.txt"
    conda: "omics"
    resources:
        cpus_per_task =5,
        threads = 5,
        runtime = 100,
        mem_mb = 20000
    shell:
        """
seqfu stats --multiqc qc/{sample}.seqfu_mqc.tsv {input.fq} > {output.stats}
        """



rule deeptools_coverage:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.filtered.bam",
    conda:
        "deeptools"
    output:
        coverage = "intermediates/{sample}/alignments/{sample}.coverage.bw"
    shell: "bamCoverage -b {input.sortedbam} -o {output.coverage} --outFileFormat bigwig"

rule deeptools_plotcoverage:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.filtered.bam",
    conda:
        "deeptools"
    output:
        genomecoverageplot = "results/{sample}/coverageGenome.png",
    shell:
        """
plotCoverage -b {input.sortedbam} \
--plotFile {output.genomecoverageplot} \
-n 1000000 \
--plotTitle "Whole Genome Coverage" \
--ignoreDuplicates \
--minMappingQuality 10
        """

rule callCoverageAndQC:
    input:
        [expand("intermediates/{sample}/alignments/coverage.bw",sample =  samples),
        expand("results/{sample}/coverageGenome.png", sample = samples),
        expand("qc/{sample}/{sample}.bamstats.txt", sample = samples)]
