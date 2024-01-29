
rule methylartistlocusplot5UTR:
    input:
        sortedbam = expand("intermediates/{sample}/alignments/{sample}.filtered.bam", sample = samples),
        bai = expand("intermediates/{sample}/alignments/{sample}.filtered.bam.bai", sample = samples)
    resources:
        runtime = 180,
        mem_mb = 5000
    params:
        ref = config["hs1sorted"],
        genes = config["genes"],
        l1hsintact = config["l1hsintact"],
        outputprefix = "results/plots/methylartist/locus/l1hsintact_"

    output:
        "outfiles/methylartistlocusplot5utr.txt"
    log:
        "logs/methylartistlocusplot5utr.log"
    conda:
        "methylartist"
    shell:
        r"""
bams=$(echo {input.sortedbam}) 2> {log}
commabams=$(echo $bams | tr ' ' ',') 2>> {log}
coords=$(awk '{{ if ($6 == "+") {{print $1":"$2"-"$2+909}} else {{print $1":"$3-909"-"$3}} }}' {params.l1hsintact})
mapfile -t coords_array <<< "$coords"
for coord in "${{coords_array[@]}}"
do
coordfnsafe=$(echo $coord | sed 's/:/-/')
if /users/mkelsey/data/tools/methylartist/methylartist locus --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt -o {params.outputprefix}${{coordfnsafe}}.png -i $coord --ref {params.ref} --motif CG -m m --highlight_bed results/tables/dmrs.bed -g /users/mkelsey/data/ref/genomes/hs1/annotations3/RTE/l1hsintact.gff.gz  -p 1,6,1,3,4 --labelgenes >> {log} 2>> {log}; then
echo "worked"
else
echo "failed"
fi
if /users/mkelsey/data/tools/methylartist/methylartist locus --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt -o {params.outputprefix}${{coordfnsafe}}.png -i $coord --ref {params.ref} --motif CG -m m --highlight_bed results/tables/dmrs.bed -g /users/mkelsey/data/ref/genomes/hs1/annotations3/RTE/l1hsintact.gff.gz  -p 1,6,1,3,4 --labelgenes --svg >> {log} 2>> {log}; then
echo "worked"
else
echo "failed"
fi
done
touch {output}       
        """

rule methylartistlocusplotexpandedview:
    input:
        sortedbam = expand("intermediates/{sample}/alignments/{sample}.filtered.bam", sample = samples),
        bai = expand("intermediates/{sample}/alignments/{sample}.filtered.bam.bai", sample = samples)
    resources:
        runtime = 180,
        mem_mb = 5000
    params:
        ref = config["hs1sorted"],
        genes = config["genes"],
        l1hsintact = config["l1hsintact"],
        outputprefix = "results/plots/methylartist/locus/l1hsintact_"

    output:
        "outfiles/methylartistlocusplotexpandedview.txt"
    log:
        "logs/methylartistlocusplotexpandedview.log"
    conda:
        "methylartist"
    shell:
        r"""
bams=$(echo {input.sortedbam}) 2> {log}
commabams=$(echo $bams | tr ' ' ',') 2>> {log}
coords=$(awk '{{print $1":"$2-6000"-"$3+6000}}' {params.l1hsintact})
mapfile -t coords_array <<< "$coords"
for coord in "${{coords_array[@]}}"
do
coordfnsafe=$(echo $coord | sed 's/:/-/')
if /users/mkelsey/data/tools/methylartist/methylartist locus --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt -o {params.outputprefix}${{coordfnsafe}}.png -i $coord --ref {params.ref} --motif CG -m m --highlight_bed results/tables/dmrs.bed -g /users/mkelsey/data/ref/genomes/hs1/annotations3/RTE/l1hsintact.gff.gz  -p 1,6,1,3,4 --labelgenes >> {log} 2>> {log}; then
echo "worked"
else
echo "failed"
fi
if /users/mkelsey/data/tools/methylartist/methylartist locus --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt -o {params.outputprefix}${{coordfnsafe}}.png -i $coord --ref {params.ref} --motif CG -m m --highlight_bed results/tables/dmrs.bed -g /users/mkelsey/data/ref/genomes/hs1/annotations3/RTE/l1hsintact.gff.gz  -p 1,6,1,3,4 --labelgenes --svg >> {log} 2>> {log}; then
echo "worked"
else
echo "failed"
fi
done
touch {output}       
        """

rule methylartistlocusplotexpandedviewltr5:
    input:
        sortedbam = expand("intermediates/{sample}/alignments/{sample}.filtered.bam", sample = samples),
        bai = expand("intermediates/{sample}/alignments/{sample}.filtered.bam.bai", sample = samples)
    resources:
        runtime = 180,
        mem_mb = 5000
    params:
        ref = config["hs1sorted"],
        genes = config["genes"],
        l1hsintact = config["l1hsintact"],
        outputprefix = "results/plots/methylartist/locus/testlocus_l1hsintact_"

    output:
        "outfiles/methylartistlocusplotexpandedviewltr.txt"
    log:
        "logs/methylartistlocusplotexpandedview.log"
    conda:
        "methylartist"
    shell:
        r"""
bams=$(echo {input.sortedbam}) 2> {log}
commabams=$(echo $bams | tr ' ' ',') 2>> {log}
coords=$(awk '{{print $1":"$2-6000"-"$3+6000}}' /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations3/RTE/ltr5intact.bed)
mapfile -t coords_array <<< "$coords"
for coord in "${{coords_array[@]}}"
do
coordfnsafe=$(echo $coord | sed 's/:/-/')
if /users/mkelsey/data/tools/methylartist/methylartist locus --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt -o {params.outputprefix}${{coordfnsafe}}.png -i $coord --ref {params.ref} --motif CG -m m --highlight_bed results/tables/dmrs.methylartistHighlight.txt -g /oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations3/RTE/ltr5intact.gff.gz  -p 1,6,1,3,4 --labelgenes >> {log} 2>> {log}; then
echo "worked"
else
echo "failed"
fi
done
touch {output}       
        """

#awk script to get gff in proper format for bgzip, tabix, and methylartist (note the space at the end and the feature being transcript)
# awk -v OFS="\t" '{print $1,"hub_3671779_hs1_hub_3671779_ncbiRefSeqCurated","transcript", $2, $3, "0.000000", $6, ".", "gene_id \"l1hsintact\"; transcript_id \"l1hsintact\"; "}' l1hsintact.bed > good.gff
# /users/mkelsey/data/tools/methylartist/methylartist locus --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt -o results/plots/methylartist/locus/testlocus_l1hsintact_chr1-71385136-71403481.png -i chr1:71384136-71404481 --ref /users/mkelsey/data/ref/genomes/hs1/hs1.sorted.fa --motif CG -m m --highlight_bed results/tables/dmrs.methylartistHighlight.txt -g /users/mkelsey/data/ref/genomes/hs1/annotations3/RTE/good.gff.gz -p 1,6,1,3,4
#awk -v OFS="\t" '{print $1,"hub_3671779_hs1_hub_3671779_ncbiRefSeqCurated","transcript", $2, $3, "0.000000", $6, ".", "gene_id \""$4"\"; transcript_id \""$4"\"; "}' ltr5intact.bed > ltr5intact.gff
rule methylartistcomposite:
    input:
        sortedbams = expand("intermediates/{sample}/alignments/{sample}.filtered.bam", sample = samples),
        bai = expand("intermediates/{sample}/alignments/{sample}.filtered.bam.bai", sample = samples)
    params:
        ref = config["hs1sorted"],
        genes = config["genes"],
        intactl1hss = config["l1hsintact"],
        l13 = config["l13"],
        l13orfs = config["l13orfs"],
        repeats = config["repeatsbgzip"]
    resources:
        runtime = 30,
        mem_mb = 100000,
        cpus_per_task = 30
    output:
       "results/plots/methylartist/composite_intactl1hs.png" 
    log:
        "logs/methylartistcomposited.log"
    conda:
        "methylartist"
    shell:
        r"""
cut -f1,2,3,6 {params.intactl1hss} > tempIntactL1hs.bed
/users/mkelsey/data/tools/methylartist/methylartist composite --bams /users/mkelsey/data/Nanopore/alz/conf/private/methylartist_bam_config.txt --meanplot_cutoff 5 -s tempIntactL1hs.bed -r {params.ref} -o {output} -t {params.l13} --mod m --primary_only -p {resources.cpus_per_task} --blocks {params.l13orfs}
        """

rule callMethylartistPlots:
    input:
        composite = "results/plots/methylartist/composite_intactl1hs.png",
        locus5utr = "outfiles/methylartistlocusplot5UTR.txt",
        locus = "outfiles/methylartistlocusplot.txt",
        locusextended = "outfiles/methylartistlocusplotexpandedview.txt"

