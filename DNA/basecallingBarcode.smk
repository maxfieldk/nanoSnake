rule dorado:
    input:
        dir = "rawdata/"
    output:
        calls = "intermediates/calls.bam"
    params:
        reference = config["reference"]
    resources:
        cpus_per_task =8,
        threads = 8,
        runtime = 5760,
        slurm_partition="gpu-he",
        mem_mb = 128000,
        slurm_extra="--time=96:00:00 --constraint=a6000 --gres=gpu:2 --mail-type=ALL --mail-user=maxfield_kelsey@brown.edu"
    shell:
        """
{params.dorado} \
basecaller \
{params.basecallingModel} \
{input.dir} \
--recursive \
--verbose \
--reference {params.reference} \
--kit-name {params.sequencingKit} > {output.calls}
        """

rule dorado_demux:
    input:
        calls = "intermediates/calls.bam"
    output:
        outfile = "outfiles/demux.txt"
    resources:
        cpus_per_task =8,
        threads = 8,
        runtime = 5760,
        mem_mb = 64000
    shell:
        """
/users/mkelsey/data/tools/dorado-0.4.3-linux-x64/bin/dorado demux --output-dir ./intermediates --no-classify {input.calls}
touch {output.outfile}
        """






