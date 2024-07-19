library(sesame)

qcs <- openSesame(snakemake@params[['idat_dir']], 
                  prep="", 
                  func=sesameQC_calcStats, 
                  BPPARAM = BiocParallel::MulticoreParam(snakemake@threads))

saveRDS(qcs, snakemake@output[[1]])
