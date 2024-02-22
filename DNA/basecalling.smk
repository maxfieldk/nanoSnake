rule dorado:
    input:
        dir = "rawdata/{rate}/{sample}"
    output:
        calls = "intermediates/{sample}/alignments/{rate}/{sample}.{type}.{modifications_string}.bam"
    wildcard_constraints:
        sample="[0-9A-Za-z]+",
        type = "[0-9A-Za-z]+",
        rate = "[0-9A-Za-z]+",
        modifications_string = "[0-9A-Za-z_-]+"
    params:
        dorado = config["dorado"],
        basecallingModel = lambda w: config["basecallingModel"][w.rate][w.type],
        reference = config["reference"]
    resources:
        cpus_per_task =12,
        threads = 12,
        slurm_partition="gpu-he",
        mem_mb = 128000,
        slurm_extra="--time=96:00:00 --constraint=a6000 --gres=gpu:2"
    shell:
        """
mkdir -p $(dirname {output.calls})
mod_string=$(echo {wildcards.modifications_string} | tr "-" ",")

{params.dorado} \
basecaller \
{wildcards.type},$mod_string \
{input.dir} \
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
