rule clair3:
    input:
        bam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam.bai"
    output:
        outfile = "outfiles/clair3AF{sample}.txt"
    params:
        ref = config["ref"],
        outdir = "intermediates/{sample}/clair3AF"
    threads: 32
    resources:
        runtime = 1000,
        mem_mb = 128000,
    conda: "clair3"
    shell:
        """
wd=$(pwd)
mkdir -p outfiles
mkdir -p {params.outdir}
run_clair3.sh \
--indel_min_af=0.05 \
--bam_fn="{input.bam}" \
--ref_fn="{params.ref}" \
--threads=24 \
--platform="ont" \
--model_path="/oscar/data/jsedivy/mkelsey/tools/remoraModels/r1041_e82_400bps_hac_v420" --output=${{wd}}/{params.outdir}
touch {output.outfile}
        """
rule clair3ROI:
    input:
        bam = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam",
        bai = "intermediates/{sample}/alignments/{sample}.sorted.filtered.bam.bai"

    output:
        outfile = "outfiles/clair3ROIAF{sample}.txt"
    params:
        ref = config["ref"],
        outdir = "intermediates/{sample}/clair3ROIAF",
        roi = "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations4/hs1.repeatMasker.l1hs.bed"
    threads: 32
    resources:
        runtime = 1000,
        mem_mb = 128000,
    conda: "clair3"
    shell:
        """
wd=$(pwd)
mkdir -p outfiles
mkdir -p {params.outdir}
run_clair3.sh \
--indel_min_af=0.05 \
--bed_fn="{params.roi}" \
--bam_fn="{input.bam}" \
--ref_fn="{params.ref}" \
--threads=24 \
--platform="ont" \
--model_path="/oscar/data/jsedivy/mkelsey/tools/remoraModels/r1041_e82_400bps_hac_v420" --output=${{wd}}/{params.outdir}
touch {output.outfile}
        """


rule sniffles:
    input:
        bam = "intermediates/{sample}/alignments/{sample}.sorted.bam"
    output:
        vcf = "intermediates/{sample}/sniffles/sniffles.vcf"
    params:
        ref = config["ref"],
        outdir = "intermediates/{sample}/sniffles"
    threads: 32
    resources:
        runtime = 1000,
        mem_mb = 128000,
    conda: "sniffles"
    shell:
        """
mkdir -p $(dirname {output.vcf})
sniffles -i {input.bam} -v {output.vcf} --threads 28 --mosaic --reference {params.ref}
        """
