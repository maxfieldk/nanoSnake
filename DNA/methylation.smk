
rule methylartistprepdssdata:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam.bai",
    params:
        ref = config["hs1sorted"],
        refindex = config["hs1sortedindex"]
    output:
        "intermediates/{sample}/{sample}.dss.tsv"
    conda: "methylartist"
    resources:
        cpus_per_task = 28,
        runtime = 1000,
        mem_mb = 80000
    shell:
        r"""
/users/mkelsey/data/tools/methylartist/methylartist wgmeth -b {input.sortedbam} --dss --ref {params.ref} --fai {params.refindex} -o {output} -p {resources.cpus_per_task} --motif CG --mod m --primary_only
        """
#be sure to put the input of an expand call in quotes when wanting it to show up properly in R script e.g. '{input.data}'
#NOTE THERE are hard coded paths and sample metadata
rule dss:
    input:
        data = expand("intermediates/{sample}/{sample}.dss.tsv", sample = samples),
    output:
        dmls = "results/tables/dmls.tsv",
        dmrs = "results/tables/dmrs.tsv"
    conda:
        "dss"
    resources:
        cpus_per_task = 10,
        mem_mb = 200000
    log:
        "logs/dss.log"
    shell:
        "Rscript scripts/dss.r data='{input.data}' dmlspath={output.dmls} dmrspath={output.dmrs} > {log} 2> test.err"

rule l1readanalysis:
    input:
        bams = expand("intermediates/{sample}/alignments/{sample}.sorted.filtered.bam", sample = samples)
    params:
        samples = samples
    output:
        l1reads = 'results/tables/l1reads.tsv'
    resources:
        cpus_per_task = 10,
        mem_mb = 50000,
        time = 200
    conda:
        "modbam"
    script:
        "scripts/readanalysis.py"

rule bedmethylanalysis:
    input:
        bedmethylpaths = expand("intermediates/{sample}/methylation/{sample}_bedMethyl.bed", sample = samples),
        dmrs = "results/tables/dmrs.tsv",
        dmls = "results/tables/dmls.tsv",
        l1reads = 'results/tables/l1reads.tsv'
    output:
        "outfiles/bedmethylanalysis.txt"
    resources:
        cpus_per_task =10,
        mem_mb = 200000,
        runtime = 400
    log:
        "logs/bedmethylanalysis.log"
    conda:
        "ds"
    shell:
        """
Rscript scripts/bedmethylanalysis.r bedmethylpaths={input.bedmethylpaths} dmrspath={input.dmrs} dmlspath={input.dmls} > {log} 2> {log}
touch {output}
        """


rule modbamtobed:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam.bai"
    params:
        ref = config["hs1sorted"],
        outdir = "results/{sample}"
    log: "logs/{sample}modbam2bed.log"
    resources:
        cpus_per_task = 40,
        mem_mb = 100000
    output:
        bed = "intermediates/{sample}/methylation/{sample}_bedMethyl.bed",
    conda:
        "omics"
    shell: 
        """
/oscar/data/jsedivy/mkelsey/tools/modkit/modkit pileup {input.sortedbam} {output.bed} \
  --ref {params.ref} \
  --preset traditional \
  -t {resources.cpus_per_task} \
  2> {log}
        """



rule modbamtobedTLDRinsertions:
    input:
        tldr = "tldr/{sample}.filtered.table.txt"
    params:
        ref = config["hs1sorted"],
        indir = "tldr/{sample}.filtered"
    resources:
        cpus_per_task = 4,
        mem_mb = 30000
    output:
        outfile = "outfiles/{sample}.bedmethyltldr.txt"
    conda:
        "omics"
    shell: 
        """
bash scripts/bedmethtldr.sh {params.indir}
touch {output.outfile}
        """

rule modbamtobedTLDRinsertionsCombinedCall:
    input:
        tldr = "tldr/SEN1.filtered_PRO1.filtered.table.txt"
    params:
        ref = config["hs1sorted"],
        indir = "tldr/SEN1.filtered_PRO1.filtered"
    resources:
        cpus_per_task = 4,
        mem_mb = 30000
    output:
        outfile = "outfiles/bedmethyltldrCombinedCall.txt"
    conda:
        "omics"
    shell: 
        """
bash scripts/bedmethtldr.sh {params.indir}
touch {output.outfile}
        """


rule callmodbamtools:
    input:
        expand("results/{sample}/outfile.txt", sample = samples)

rule modbamtools:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.bam",
    params:
        refseq = config["genes"],
        outdir = "results/{sample}"
    output:
        outfile = "results/{sample}/outfile.txt"        
    conda:
        "modbam"
    shell:
        """
modbamtools plot -r chr20:58815000-58895000 \
    --gtf {params.refseq} \
    --out {params.outdir} \
    --prefix {wildcards.sample} \
    --samples {wildcards.sample} \
    --track-titles Genes \
    {input.sortedbam}
touch {output.outfile}
        """
rule extractmods:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
    params:
        ref = config["hs1sorted"],
        outdir = "results/{sample}"
    log: "logs/{sample}modbam2bed.log"
    resources:
        cpus_per_task = 16,
        mem_mb = 100000
    output:
        bed = "intermediates/{sample}/methylation/{sample}_mods.tsv",
    conda:
        "omics"
    shell: 
        """
/oscar/data/jsedivy/mkelsey/tools/modkit/modkit extract \
{input.sortedbam} \
{output.bed} \
--region chr20 \
--ref {params.ref} \
-t {resources.cpus_per_task} \
2> {log}
        """

rule callmodbamtobedandextractmods:
    input:
        [expand("intermediates/{sample}/methylation/{sample}_bedMethyl.bed", sample = samples),
        expand("intermediates/{sample}/methylation/{sample}_mods.tsv", sample = samples)]
