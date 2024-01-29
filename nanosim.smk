rule nanosimCharaceterization:
    input:
        aln = "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/rna_7/alignments/genome/rna_7.sorted.bam",
        reads = "/users/mkelsey/data/Nanopore/dRNALF1/intermediates/rna_7/fastqs/rna_7.fq.gz"
    output:
        outfile =  "outfiles/nanosim.out",
    params:
        reference = "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.nonrefcontigs.fa",
        annotation = "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations4/hs1.repMask.refseq.nonref.sorted.gtf",
        transcriptome = "/users/mkelsey/data/ref/genomes/hs1/annotations4/hs1.transcripts.repMask.refseq.nonref.fna.fai"
    resources:
        cpus_per_task =32,
        mem_mb = 128000
    conda:
        "nanosim"
    shell:
        """
mkdir -p $(dirname {output.outfile})
read_analysis.py transcriptome -rg {params.reference} -rt {params.transcriptome} -annot {params.annotation} -ga {input.aln} -i {input.reads} -o $(dirname {output.outfile}) -t 28 -c
        """
        

