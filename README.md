# methylation_array

methylation_array is a Snakemake workflow designed to handle Illumina methylation array data through mutliple steps: QC, creating SummarizedExperiment object, generating PCA plot, and cconduct differential methylation analysis. The workfow uses R bioconductor package [SeSAME](https://bioconductor.org/packages/release/bioc/html/sesame.html) as the central component. The workflow is optimized for use on the HPC system at Van Andel Insititute, but it can also be adapted for use on any Unix-based system with some adjustments to the environment setup.

## Installation

No installation is required for using this, you need to git clone this repo. 
```
git clone https://github.com/vari-bbc/methylation_array.git
```

## Configure the workflow
  * Step 1. Move your idat files to directory "idat_files". Each of your sample should have two idat files, one for green channel and one for red channel.   
  * Step 2. List the names of these idat files in bin/units.tsv. See the example comes with the repo, do not change the column name "idat". 
  * Step 3. Copy the SampleSheet.csv, which you often receive from the sequencing facility, under bin/. If this sample sheet has a different name, change it to "SampleSheet.csv". Open SampleSheet.csv, and count how many rows were header rows. The sample information often starts below "Sample_Name", count how many rows were above "Sample_Name". You will need this information below. 
  * Step 4. Fill the config.yaml file under bin/. samplesheet_header_rows is the number from step 3. array_type can be EPICv2/EPIC/HM450 for human, and MM285 for mouse. Recommonded prep_method for EPICv2/EPIC/HM450 and MM285 are QCDPB and TQCDPB respectively. For more information please refer to [SeSAMe vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html) on bioconductor. modules can be left alone. 
  * Step 5 Fill the meta.tsv file under bin/. This file provides important information to the workflow as it links the idat files to your samples. The 5th column, Sample_Name should have the same sample_names from SampleSheet.csv (1st column of SampleSheet.csv). Group column (1st column in meta.tsv) is important, this should be grouping information for you samples in differential methylation analysis. It can be the same as genotype, or diet, populations etc. Or it can be a combination of genotype and diet. For instance you have a 2 x 2 design experiment of two genotypes, WT and KO, and two diet levels, high fat diet (hfd) and standard chow diet (scd), you want to look at methylation difference between KO_hfd and KO_scd; for this case, you will put KO_hfd, KO_scd, WT_hfd and WT_scd as group informaiton. Other columns can be left empty or you can fill with the same value for all the samples.  
  * Step 6. Lastly, fill comparisons.tsv; using example from step 5, you can put KO_hfd as group1 and KO_scd as group2. group2 will be your reference group, this means when we calculate differential methylation, it is group1-group2, postive value means methylation in group1 is higher than group2. You can list more than one comparison, using the above example, you can also look at KO_hfd and KO_scd. Make sure each compaprison is in one row.  

## Run the workflow
After you the setup you have done, you are ready to run the workflow. Make sure to go back to the main dir where Snakefile and README file are located. Then run:
```
sbatch bin/run_snake.sh
```   

Or if you want to do a dry run and see what rules will be executed, you can do:
```
module load bbc2/snakemake
snakemake -npr
```


## Results
The workflow will generate QC results, Summerizedexperiment object, DML, DMP results, a PCA plot and volcano plots for each comparison. 
  * QC results is a named list, names correspond to chip and position ID of each sample
  * SummarizedExperiment object (se.rds) has beta values for all samples (it can be accessed using assay(se), or assays(se)$betas), as well as meta information(it can be accessed using colData(se))
  * Differential methylation results includes four components: SummarizedExperiment, dml_smry, dml_result, and dmr result. SummarizedExperiment could be useful for plotting heatmap. dml_smry has the information of linear model fit for each CpG site, you do not usually need this file. dml_result is the summary statisics of dml_smry, one row per CpG site. Note, the pvalue here is not adjusted for mulitple comparison. Lastly, dmr is the results for differential methylated regions (DMR). DMR merges neighboring CpG sites that show consistent methylation variation into a region. A DMR can contains multiple CpG sites or only one CpG site. 
  * Two types of plots are generated, they are under "analysis/plots". PCA is a dimensional reduction of the beta values into PC1 and PC2, the data points are color coded based on the group information from meta.tsv. For dml, a volcano plot is generated for each comparison. Each CpG site is represents as a dot in the plot. Orange and blue color indicate that CpG site passed p-value cutoff (default = 0.05). Red and blue color indicate that CpG site passed beta value difference threshold (default = 0.05 or 5%). This means at least 5% of the methylation level difference between two groups you were comparing. Grey color indicates that CpG site did not pass either of beta value difference or pvalue threshold. Note, the p-value in volcano plot is not adjusted for multiple comparison. It provides a good visualzation of the methylation trend between two groups.   
