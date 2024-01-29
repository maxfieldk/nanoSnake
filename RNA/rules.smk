####
# Good tools to make new rules out of - filtlong - seqfu - 
rule dorado:
    input:
        dir = "rawdata/{sample}"
    params:
        dorado = config["dorado"],
        basecallingModel = config["basecallingModel"],
        reference = config["reference"]
    output:
        fqcalls = "intermediates/{sample}/fastqs/{sample}.fq.gz"
    resources:
        cpus_per_task =8,
        threads = 8,
        runtime = 5760,
        slurm_partition="gpu-he",
        mem_mb = 128000,
        slurm_extra="--time=96:00:00 --gres=gpu:2 --mail-type=ALL --mail-user=maxfield_kelsey@brown.edu"
    shell:
        """
{params.dorado} \
basecaller \
{params.basecallingModel} \
{input.dir} \
--recursive \
--verbose \
--emit-fastq | gzip > {output.fqcalls}
        """
# --constraint=a6000

rule isoquant:
    input:
        fqcalls = expand("intermediates/{sample}/fastqs/{sample}.fq.gz", sample = samples)
    output:
        "intermediates/isoquantNEW/done.out"
    priority: 100
    params:
        reference = config["reference"],
        annotation = config["annotation"]
    conda: "isoquant"
    resources:
        cpus_per_task =32,
        runtime = 4000,
        mem_mb = 128000,
    shell:
        """
isoquant.py -d nanopore --stranded forward --fastq {input.fqcalls} \
 --reference {params.reference} --genedb {params.annotation} \
 --output $(dirname {output}) --threads 30
 touch {output}
        """

rule isoquantresume:
    input:
        fqcalls = expand("intermediates/{sample}/fastqs/{sample}.fq.gz", sample = samples)
    output:
        "intermediates/isoquantNEW/done.resume.out"
    priority: 100
    conda: "isoquant"
    resources:
        cpus_per_task =32,
        runtime = 4000,
        mem_mb = 128000,
    shell:
        """
isoquant.py --resume --output $(dirname {output}) --threads 30
touch {output}
        """

rule fq_to_DNA_fasta_for_repeatmasker:
    input:
        fq = "intermediates/{sample}/fastqs/{sample}.fq.gz"
    output:
        fa = "intermediates/{sample}/fastqs/{sample}.fa"
    resources:
        cpus_per_task =20,
        runtime = 900,
        mem_mb = 80000
    conda:
        "omics"
    shell:
        """
seqkit fq2fa {input.fq} > {output.fa}.temp
seqkit seq --rna2dna {output.fa}.temp > {output.fa}
rm {output.fa}.temp
        """

def inputToMergeGTFs(wildcards):
    checkpoint_output = checkpoints.split_reads_fasta.get(**wildcards).output[0]
    sample = wildcards.sample
    print("sample is %s"%sample)
    part_nums=glob_wildcards(os.path.join(checkpoint_output, "%s.part_{part_num}.fa"%sample)).part_num
    print("part_nums is %s"%part_nums)
    expand_call = expand("intermediates/{{wildcards.sample}}/fastqs/repeatmasker/{part_num}/{{wildcards.sample}}.part_{part_num}.fa.gtf",part_num = part_nums)
    print("expand_call is %s"%expand_call)
    list_call = ["intermediates/%s/repeatmasker/%s/%s.part_%s.fa.gtf"%(sample, part_num, sample, part_num) for part_num in part_nums]
    print("list_call is %s"%list_call)
    return(list_call)


checkpoint split_reads_fasta:
    input:
        fa = "intermediates/{sample}/fastqs/{sample}.fa"
    output:
        directory("intermediates/{sample}/fastqs/split")
    conda:
        "omics"
    shell:
        """
mkdir -p {output}
seqkit split2 {input.fa} -p 4 -f --out-dir {output}
        """

rule repeatmasker:
    input:
        chr_fasta = "intermediates/{sample}/fastqs/split/{sample}.part_{part_num}.fa"
    output:
        rmout = "intermediates/{sample}/repeatmasker/{part_num}/{sample}.part_{part_num}.fa.out"
    resources:
        cpus_per_task =20,
        runtime = 1200,
        mem_mb = 50000
    conda:
        "omics"
    shell:
        """
mkdir -p $(dirname {output})
/oscar/data/jsedivy/mkelsey/tools/RepeatMasker/RepeatMasker -species human -pa {resources.cpus_per_task} -gff {input.chr_fasta} -dir $(dirname {output})
        """


rule getGtfs:
    input:
        rmout = "intermediates/{sample}/repeatmasker/{part_num}/{sample}.part_{part_num}.fa.out"
    output:
        gtfout = "intermediates/{sample}/repeatmasker/{part_num}/{sample}.part_{part_num}.fa.gtf"
    conda:
        "evo"
    shell:
        """
scripts/outToGtf.sh {input.rmout} {output.gtfout}
        """

rule mergeGtfsandCleanupRM:
    input:
        inputToMergeGTFs
    output:
        gtf = "intermediates/{sample}/repeatmasker/{sample}_repeatmasker.gtf"
    conda:
        "evo"
    shell:
        """
cat {input} > {output.gtf}
sort -k1,1V -k4,4n -k5,5n {output.gtf} > tmp.gtf
ls -d RM_* | xargs rm -r
mv tmp.gtf {output.gtf}
        """


# rule repeat_mask_reads:



# RepeatMasker -pa 6 -dir “output_directory” -nolow -norna -species human -no_is -a -u -xsmall -xm “LINE_input_reads.fasta”




rule minimap2RNAGENOME:
    input:
        fq = "intermediates/{sample}/fastqs/{sample}.fq.gz"
    params:
        reference = config["reference"],
        annotation = config["annotation"],
        junctionbed = config["junctionbed"]
    output:
        bam = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam"
    resources:
        cpus_per_task =20,
        runtime = 900,
        mem_mb = 80000
    conda:
        "minimap2"
    shell:
        """
minimap2 -ax splice -uf -N 10 -k14 -t 12 --junc-bed {params.junctionbed} {params.reference} {input.fq} | \
samtools sort -@4 -T {wildcards.sample} -O bam -o {output.bam}
samtools index -@8 {output.bam}
samtools stats {output.bam} > {output.bam}.stats.txt        
        """

rule minimap2RNATRANSCRIPTOME:
    input:
        fq = "intermediates/{sample}/fastqs/{sample}.fq.gz"
    params:
        referencetranscriptome = config["referencetranscriptome"]
    output:
        bam = "intermediates/{sample}/alignments/transcriptome/{sample}.sorted.bam"
    resources:
        cpus_per_task =20,
        runtime = 900,
        mem_mb = 80000
    conda:
        "minimap2"
    shell:
        """
minimap2 -ax map-ont -uf -N 10 -k14 -t 12  {params.referencetranscriptome} {input.fq} | \
samtools sort -@4 -T {wildcards.sample} -O bam -o {output.bam}
samtools index -@8 {output.bam}
samtools stats {output.bam} > {output.bam}.stats.txt        
        """

rule NanoCountgenesandRTEs:
    input:
        sortedbam = "intermediates/{sample}/alignments/transcriptome/{sample}.sorted.bam"
    output:
        counts = "intermediates/{sample}/counts/transcriptome/{sample}.nanocount.tsv",
    conda:
        "omics"
    resources:
        cpus_per_task =4,
        runtime = 3000,
        mem_mb = 32000,
    shell: 
        """
NanoCount -i {input.sortedbam} -o {output.counts} --extra_tx_info
        """

rule featureCountsgenesandRTEs:
    input:
        sortedbam = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam"
    output:
        countsmessy = "intermediates/{sample}/counts/genome/{sample}rtesandgenes.counts_messy.txt",
        counts = "intermediates/{sample}/counts/genome/{sample}rtesandgenes.counts.txt",
    params: 
        annotation = config["annotation"]
    conda:
        "omics"
    resources:
        cpus_per_task =4,
        runtime = 3000,
        mem_mb = 32000,
    shell: 
        """
featureCounts -L -M --primary --ignoreDup --largestOverlap -R CORE -a {params.annotation} -o {output.countsmessy} {input.sortedbam}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
        """

rule featureCountsgenesandRTEsSTRINGENT:
    input:
        sortedbam = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam"
    output:
        countsmessy = "intermediates/{sample}/counts/genome/stringent/{sample}rtesandgenes.counts_messy.STRINGENT.txt",
        counts = "intermediates/{sample}/counts/genome/stringent/{sample}rtesandgenes.counts.STRINGENT.txt",
    params: 
        annotation = config["annotation"]
    conda:
        "omics"
    resources:
        cpus_per_task =4,
        runtime = 3000,
        mem_mb = 32000,
    shell: 
        """
mkdir -p intermediates/{wildcards.sample}/counts/genome/stringent
featureCounts -L -M --primary --ignoreDup --largestOverlap --fracOverlapFeature 0.8 --fracOverlap 0.8 -R CORE -a {params.annotation} -o {output.countsmessy} {input.sortedbam}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
        """



rule filterbam:
    input:
        bam = "intermediates/{sample}/alignments/{sample}.sorted.bam"
    output:
        bam = "intermediates/{sample}/alignments/{sample}.filtered.bam"
    resources:
        cpus_per_task =6,
        mem_mb = 40000
    conda:
        "omics"
    log: "logs/{sample}_filterbam.log"
    shell: "samtools view -b -F 0x100 -q 10 -e '[qs] > 10' {input.bam} > {output.bam} 2> {log}"

rule sortbam:
    input:
        bam = "intermediates/{sample}/alignments/{directory}/{sample}.bam"
    output:
        sortedbam = "intermediates/{sample}/alignments/{directory}/{sample}.sorted.bam"
    resources:
        cpus_per_task =10,
        mem_mb = 128000
    conda:
        "omics"
    shell: "samtools sort -@8 -m4g {input.bam} > {output.sortedbam} 2> {log}"
    
rule indexbam:
    input:
        sortedbam = "intermediates/{sample}/alignments/{directory}/{bamname}.bam"
    output:
        bai = "intermediates/{sample}/alignments/{directory}/{bamname}.bam.bai"
    resources:
        cpus_per_task =10,
        mem_mb = 128000
    conda:
        "minimap2"
    shell: "samtools index  -@6 {input.sortedbam}"

rule bamstats:
    input:
        sortedbam = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam",
        bai = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam.bai"
    output:
        "qc/{sample}/{sample}.bamstats.txt"
    conda:
        "omics"
    shell:
        """
samtools stats {input.sortedbam} > {output}
        """

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
    
rule dorado_seqsummary:
    input:
        sortedbam = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam",
        bai = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam.bai"
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
        sortedbam = "intermediates/{sample}/alignments/genome/{sample}.sorted.bam",

    output:
        "qc/{sample}/{sample}pycoQC.html"
    conda:
        "pycoQC"
    shell:
        """
pycoQC --summary_file {input.seqsummary} --bam_file {input.sortedbam} --html_outfile {output} --min_pass_qual 10 --sample
        """


rule multiqc:
    conda:
        "omics"
    params:
        outputdir = "qc"
    output:
        report = "qc/multiqc_report.html",
    shell:
        """
multiqc --force --outdir {params.outputdir} --filename {output.report} --export .
        """