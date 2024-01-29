

rule dorado:
    input:
        dir = "rawdata_4khz/{sample}"
    output:
        calls = "intermediates/{sample}/alignments/{sample}.bam"
    params:
        dorado = config["dorado"],
        basecallingModel = lambda w: config["basecallingModel"]["4khz"],
        reference = config["reference"]
    resources:
        cpus_per_task =12,
        threads = 12,
        slurm_partition="gpu-he",
        mem_mb = 128000,
        slurm_extra="--time=96:00:00 --constraint=a6000 --gres=gpu:2"
    shell:
        """
{params.dorado} \
basecaller \
{params.basecallingModel} \
{input.dir} \
--modified-bases  5mCG_5hmCG \
--recursive \
--verbose \
--reference {params.reference} > {output.calls}
        """



    # resources:
    #     cpus_per_task =8,
    #     threads = 8,
    #     runtime = 5760,
    #     slurm_partition="gpu-he",
    #     mem_mb = 128000,
    #     slurm_extra="--time=96:00:00 --constraint=a6000 --gres=gpu:2 --mail-type=ALL --mail-user=maxfield_kelsey@brown.edu"
