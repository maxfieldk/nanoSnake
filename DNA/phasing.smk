
rule whatshap:
    input:
        shortreadsbam = "../../freebayes/alignment/lf1.sorted.bam",
        shortreadsbamindex = "../../freebayes/alignment/lf1.sorted.bam.bai",
        longreadsbam = "intermediates/{sample}/alignments/{sample}.sorted.bam",
        longreadsbamindex = "intermediates/{sample}/alignments/{sample}.sorted.bam.bai",
        vcf = "../../freebayes/results/variants.vcf",
        ref = config["hs1sorted"],
    output:
        phasedvcf = "intermediates/{sample}/alignments/phased.vcf"
    resources:
        runtime = 1200
    log:
        "logs/whatshap{sample}.log"
    conda:
        "whatshap"
    resources:
        cpus_per_task =10
    shell:	
        """
whatshap phase \
--ignore-read-groups \
-o {output.phasedvcf} \
--reference={input.ref} \
{input.vcf} \
{input.shortreadsbam} \
{input.longreadsbam} > {log}
        """

rule bgzipandindexphased:
    input:
        phasedvcf = "results/phased.vcf"
    output:
        phasedvcfgz = "results/phased.vcf.gz"
    log:
        "logs/bgzipandindexphased.log"
    conda:
        "omics"
    resources:
        cpus_per_task =10
    shell:	
        """
bgzip {input.phasedvcf}
tabix {input.phasedvcf}.gz
        """
        
rule addphasetobam:
    input:
        sortedbam = "intermediates/{sample}/alignments/{sample}.sorted.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.bam.bai",
        phasedvcfgz = "../../freebayes/results/phased.vcf.gz",
        ref = config["hs1sorted"]
    log:
        "logs/addphasetobam{sample}.log"
    output:
        haplotaggedbam = "intermediates/{sample}/alignments/{sample}.sorted.haplotagged.bam",
    resources:
        cpus_per_task = 20
    conda:
        "whatshap"
    shell:
        """
whatshap haplotag -o {output.haplotaggedbam} --output-threads=$(( {resources.cpus_per_task}-10 )) --ignore-read-groups --reference {input.ref} {input.phasedvcfgz} {input.sortedbam} 2> {log}
        """