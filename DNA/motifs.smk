rule gimme:
    input:
        hypo = "Rintermediates/promoters_dmhyporegions.bed",
        hyper = "Rintermediates/promoters_dmhyperregions.bed"
    output:
        "outfiles/gimme.txt"
    resources:
        cpus_per_task = 26,
        mem_mb = 64000,
        runtime = 300
    log:
        "logs/gimme.log"
    conda:
        "Rclusterprofiler"
    shell:
        """
bedtools getfasta -fi /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa -bed Rintermediates/promoters.bed -fo Rintermediates/promoters.fa
bedtools getfasta -fi /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa -bed Rintermediates/promoters_dmhyperregions.bed -fo Rintermediates/promoters_dmhyperregions.fa
bedtools getfasta -fi /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa -bed Rintermediates/promoters_dmhyporegions.bed -fo Rintermediates/promoters_dmhyporegions.fa
gimme motifs results/tables/dmrs_hyper.bed results/gimme/hyper --background gc -g /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa --known --nthreads 24 -p /oscar/data/jsedivy/mkelsey/ref/pwms/hocomoco_jasperformat.txt --size 500 --noreport
gimme motifs results/tables/dmrs_hypo.bed results/gimme/hypo --background gc -g /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/hs1.sorted.fa --known --nthreads 24 -p /oscar/data/jsedivy/mkelsey/ref/pwms/hocomoco_jasperformat.txt --size 500 --noreport
touch {output}
        """
        
rule homer_find_motifs:
    input:
        dmrs = "results/tables/dmrs.tsv"
    params:
        hs1 = config["hs1sorted"],
        outputdir = "results/homer"
    output:
        "homeroutfile.txt"
    conda:
        "omics"
    log:
        "logs/homer_find_motifs.log"
    shell:
        r"""
awk '{{IFS=OFS="\t"}}; {{print "DMR" NR, $1, $2, $3, "+"}}' {input.dmrs} > tempbed.bed
findMotifsGenome.pl tempbed.bed {params.hs1} {params.outputdir} -size given
touch {output}
        """