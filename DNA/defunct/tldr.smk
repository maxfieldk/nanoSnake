rule tldr:
    input:
        bam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam.bai"
    output:
        outfile = "outfiles/tldr.{sample}.outfile.txt"
    resources:
        cpus_per_task = 28,
        mem_mb = 100000,
        runtime = 300
    conda:
        "tldr"
    shell:
        """
mkdir -p tldr
cd tldr
tldr -b ../{input.bam} \
-e /oscar/data/jsedivy/mkelsey/tools/tldr/ref/teref.ont.human.fa \
-r /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa \
-p 24 \
--detail_output \
--extend_consensus 4000 \
--trdcol \
--color_consensus \
--flanksize 500 \
--keep_pickles \
--methylartist
cd ..
touch {output.outfile}
        """


rule tldrALLSAMPLES:
    input:
        bam = expand("intermediates/{sample}/alignments/{sample}.sorted.filtered.bam", sample = samples),
        bai = expand("intermediates/{sample}/alignments/{sample}.sorted.filtered.bam.bai", sample = samples)
    output:
        outfile = "outfiles/tldrall.outfile.txt"
    resources:
        cpus_per_task = 28,
        mem_mb = 100000,
        runtime = 300

    conda:
        "tldr"
    shell:
        """
mkdir -p tldr
cd tldr
bams1=({input.bam})
bams2=$(printf "\"../%s\"," "${{bams1[@]}}")
bams=${{bams2%?}}
tldr -b $bams \
-e /oscar/data/jsedivy/mkelsey/tools/tldr/ref/teref.ont.human.fa \
-r /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa \
-p 24 \
--detail_output \
--extend_consensus 500 \
--trdcol \
--color_consensus \
--flanksize 500 \
--keep_pickles \
--methylartist
cd ..
touch {output.outfile}
        """