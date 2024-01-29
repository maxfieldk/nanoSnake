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
    