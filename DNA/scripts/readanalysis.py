import logging
import pysam
import pandas as pd
import modbampy
from modbampy import ModBam
import matplotlib as mpl
import numpy as np
import yaml

print("script is running")
with open('conf/config.yaml', 'r') as file:
    conf = yaml.safe_load(file)

intactl1srefnonrefpath = conf["intactl1srefnonref"]
l1hsintactpath = conf["l1hsintact"]
l1hspid80path =  conf["l1hspid80"]
l1pa2pid80path = conf["l1pa2pid80"]
sample_tablepath = conf["sample_table"]

sample_table = pd.read_csv(sample_tablepath)

try:
    bampaths = snakemake.input["bams"]
    l1readspath = snakemake.output["l1reads"]
    samples = snakemake.params["samples"]
    print("snakemake variables loaded")
except:
    print("snakemake variables not loaded")
    samples = sample_table.sample_name
    l1readspath = 'results/tables/l1reads.tsv'
    bampaths = ["intermediates/%s/alignments/%s.sorted.filtered.bam"%(sample, sample) for sample in samples]

intactl1srefnonref = pd.read_csv(intactl1srefnonrefpath, delimiter="\t")
l1hsintact = pd.read_csv(l1hsintactpath, delimiter="\t", header=None)
l1hspid80 = pd.read_csv(l1hspid80path, delimiter="\t", header=None)
l1pa2pid80 = pd.read_csv(l1pa2pid80path, delimiter="\t", header=None)


l1hsintact["type"] = "l1hsintact"
l1hspid80["type"] = "l1hspid80"
l1pa2pid80["type"] = "l1pa2pid80"
l1 = pd.concat([l1hsintact, l1hspid80, l1pa2pid80])

l1 = intactl1srefnonref
##############

l1dict = {}
for sample in samples:
    l1dict[sample] = {}
    condition = sample_table.loc[sample_table["sample_name"] == sample, "condition"].iloc[0]
    for index, row in l1.iterrows():
        print(row["chr"])
        reads = []
        with ModBam("intermediates/%s/alignments/%s.sorted.filtered.bam"%(sample, sample)) as bam:
            for read in bam.reads(row["chr"], row["start"], row["end"]):
                for pos_mod in read.mod_sites:
                    ele = [*pos_mod]
                    reads.append(ele)
        try:
            df = pd.DataFrame(reads, columns=['id', 'refpos',
                                            'querypos', 
                                            'modstrand', 'unclear',
                                            'canon_base',
                                            'mod_base', 'mod_base_score'])
            df["element_strand"] = row["strand"]
            df["element_start"] = row["start"]
            df["element_stop"] = row["stop"]
            df["element_chr"] = row["chr"]
            df["element_type"] = row["family"]
            df["element_refstatus"] = row["source"]
            df["element_intactness"] = row["intactness"]
            df["element_uid"] = "{}_{}_{}_{}".format(df["element_chr"][0], df["element_start"][0], df["element_stop"][0], df["element_strand"][0])
            df = df.assign(meth = np.where(df.mod_base_score > 255/2, 1, 0))
            df["condition"] = condition
            df["sample"] = sample
            l1dict[sample][df["element_uid"][0] + "_" + df["element_type"][0]] = df
        except:
            print("exception", next(read.mod_sites))

valuesextracted = [pd.concat(l1dict[sample].values(), axis=0) for sample in samples]
concatenated_dataframe = pd.concat(valuesextracted, axis=0)

concatenated_dataframe.to_csv(l1readspath, sep='\t', index=False)