import pandas as pd
import os
from shutil import which
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.0.0")

configfile: "bin/config.yaml"

#units = pd.read_table("bin/units.tsv")
comps = pd.read_table("bin/comparisons.tsv")

##### target rules #####

rule all:
    input:
        "analysis/sesame_se/se.rds",
        "analysis/sesame_qc/qc.rds",
        expand("analysis/sesame_DM/{comp.group1}.v.{comp.group2}.rds", comp=comps.itertuples()),

rule sesame_qc:
    input:
        idat_dir=config['idat_dir']
        #expand("idat_files/{file}", file=units.idat) 
    output:
        "analysis/sesame_qc/qc.rds"
    benchmark:
        "benchmarks/sesame_qc/bench.txt"
    #params:
     #   idat_dir=config['idat_dir']
    threads: 8
    resources:
        mem_gb=96,
        log_prefix="sesame_qc" 
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/run_sesame_QC.R"

rule sesame_se:
    input:
        idat_dir=config['idat_dir'],
        #idats = expand("idat_files/{file}", file=units.idat),
        samplesheet = "bin/SampleSheet.csv",
        meta = "bin/meta.tsv"
    output:
        "analysis/sesame_se/se.rds"
    benchmark:
        "benchmarks/sesame_se/bench.txt"
    params:
       # idat_dir="idat_files/",
        manifest_id=config['sesame_params']['array_type'],
        sesame_prep=config['sesame_params']['prep_method'], # recommended by vignette for MM285
        samplesheet_skip=config['samplesheet_header_rows']
    threads: 8
    resources:
        mem_gb=96,
        log_prefix="sesame_se"
    envmodules:
        config['modules']['R']
    script: 
        "bin/scripts/make_sesame_SE.R"


rule sesame_DM:
    input:
        "analysis/sesame_se/se.rds"
    output:
        "analysis/sesame_DM/{group1}.v.{group2}.rds"
    benchmark:
        "benchmarks/sesame_DM/{group1}.v.{group2}.txt"
    params:
        array_type=config['sesame_params']['array_type']
    threads: 4
    resources:
        mem_gb=72,
        log_prefix=lambda wildcards: "_".join(wildcards)
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/run_sesame_DM.R"


