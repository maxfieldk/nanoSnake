rule subsetFqByReadIDsForAssembly:
    input:
        fq = "intermediates/{sample}/fastqs/{sample}.fq",
        readIDs = "intermediates/{sample}/alignments/{sample}_to_{target}.readIDs.txt"
    output:
        fq = "results/assembly/{sample}/forAssembly_{sample}_to_{target}forAssembly_.fq",
    conda: "omics"
    resources:
        cpus_per_task =10,
        threads = 10,
        runtime = 100,
        mem_mb = 32000,
    shell:
        """
mkdir -p $(dirname {output.fq})
awk '!seen[$0]++' {input.readIDs} > {input.readIDs}.duplicatefiltered
seqtk subseq {input.fq} {input.readIDs}.duplicatefiltered  > {output.fq}
rm {input.readIDs}.duplicatefiltered
        """

rule flyeAll:
    input:
        fq = "/users/mkelsey/data/Nanopore/l1kd/intermediates/L1KD/fastqs/L1KD.fq"
    output:
        outfile =  "assembly/all.out",
    resources:
        cpus_per_task =32,
        mem_mb = 200000,
        runtime = 2400
    conda:
        "omics"
    shell:
        """
mkdir -p $(dirname {output.outfile})
flye --nano-raw {input.fq} --out-dir $(dirname {output.outfile}) --threads 30
        """

rule flye:
    input:
        fq = "results/assembly/{sample}/forAssembly_{sample}_to_{target}forAssembly_.fq"
    output:
        outfile =  "results/assembly/{sample}/{target}/{target}.out",
    resources:
        cpus_per_task =32,
        mem_mb = 128000
    conda:
        "omics"
    shell:
        """
mkdir -p $(dirname {output.outfile})
flye --nano-raw {input.fq} --out-dir $(dirname {output.outfile}) --threads 30
        """

