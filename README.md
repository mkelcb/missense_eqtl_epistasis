# missense X eQTL epistasis detection


Pipeline overiview of 3 stages:

**I) Extract raw data for the missense/eQTL SNPs and phase them into .haps format**
1) This involves preparing data for the annotations from the de Lange/Huang studies' supplementary tables, and feeding them though a few bash/python scripts and SHAPEIT3 
2) I have not included these, as these will be completely rewritten to only include fine mapped missense SNPs and also try to get more eQTL SNPs with lower p-values.
3) For illustration purposes, I just added the final .haps files (data/SNPs/rs2476601/) that can be used in Stages I and II. (they were phased together with ~1000 SNP supporting variants flanking them)

**II) Produce design matrices for the models:**
1) *python scripts/produceDesignMatrix-Chris.py --annotationFileLoc data/misc/missense_eqtl_annotation.txt --baseDir data/SNPs/*
- this will output the two design matrices for the #BAD_HAPLO and the haplo based regression models (that is yours)
- the code for this is very straightforward python, should be easy to understand
- this also includes how the 'BAD_HAPLO' predictor is generated

**III) Model fit:**
1) *Rscript scripts/modelFit-Chris.R data/SNPs/ data/misc/covarOut_covsAll 0.005 results/ data/misc/missense_eqtl_annotation.txt* data/misc/gwas1_3_ibdseq_ukbb_ukbb_merged_all_subtypes.phe 0 data/misc/missense_gene_lookup.csv 0 0
- this fits the model in different ways (depending on arguments): phenotype specific (IBD, UC, CD or UC vs CD) and it produces results for each different cohort. Each of the 5 cohorts are fit separately, then combined, and also a meta-analysed via Stouffer's method. It can also produce marginal models.
- While I tried to make this script as readable as possible as well, there are many conditional sections depending on the phenotype, cohorts and marginal models that are repetitive but necessary. Also the bits where I match the covariates to the haplotype specific regression (your) model is also quite tedious.
2) I saved an R session under results/workspace.RData that has most objects loaded in
3) As you can see, I already fit your regression model too, I just haven't applied the 'post-processing' step of fixing the inflated variances, as I just never had the time to implement that. So if you want to add that, the 2 places where that would go can be found by searching for the string *"### ADD VARIANCE CORRECTION TO THIS MODEL HERE ###"*
4) Once we change the datasets the 'Combined' analysis will be removed, and instead we will do a fixed-effects meta analysis
5) I did not include the code for producing tables for the Genotype/Haplotype frequencies. If you want that too let me know.
