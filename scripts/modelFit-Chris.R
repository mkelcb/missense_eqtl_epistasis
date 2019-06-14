#!/usr/bin/env Rscript

# Fits regression models onto a design matrices produced by produceDesignMatrix-Chris.py
# example command:
# Rscript scripts/modelFit-Chris.R data/SNPs/ data/misc/covarOut_covsAll 0.005 results/ data/misc/missense_eqtl_annotation.txt data/misc/gwas1_3_ibdseq_ukbb_ukbb_merged_all_subtypes.phe 0 data/misc/missense_gene_lookup.csv 0 0


require("metap")
if (!require("metap")) install.packages("metap")
library("metap")
if (!require("flextable")) install.packages("flextable")
library("flextable")

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 10) {stop("not enough Arguments received")} 

# Process command line arguments:
baseDir= args[1] # where to start looking for subfolders
covsFileLoc = args[2]  # the covariates file in plink format with a header, and one row for each 
sigThreshold = args[3] # the significance threshold to be used 
outLoc = args[4] # where the output should be saved
annotationDataLoc = args[5]
uc_vs_cdLoc = args[6] # Controls = 1 = UC, Cases = 2 = CD
subPhenoMode =args[7] # if we only look at UC (1) CD (2) or vs each other (3), or if all (0)
geneNameLookupLoc = args[8] # a separate annotation file to be able to link the missenseRS ids to a gene symbol (Carl asked for this)
otherPCsAddedAsCovars = (args[9] == 1) # if PCs should be used other than PC1
MARGINAL= (args[10] == 1) # if marginal/univariate regression should be fit instead of multiple regression


##############################################################
# Helper functions
##############################################################

# extracts p-values from a glm model, even when there are NA terms, it returns the p-values in a vector, with the correct index for each temr
getPvalues = function (model) {   #  model = myModel_glm_current
  coeficients = coef(model, complete = FALSE)
  numTerms = length(coeficients)
  extractedPvals = summary(model)$coefficients[,4]
  
  # if there are less P-values than terms, then we had an NA somewhere
  if (length( extractedPvals) != numTerms) {
    NA_indices = which(is.na(coeficients)) # get the indices where we have NA
    for( k in 1:length(NA_indices) ) { 
      extractedPvals = append(extractedPvals, NA, after = NA_indices[k] -1)
    }
  }
    return(extractedPvals)
}


# fits a univariate model for a design matrix, where each term is fit separately from the others to obtain marginal p-values
getMarginalPvals = function (current_cohort) {  
  y = current_cohort$y
  pvals = c()
  for( k in 2:ncol(current_cohort) ) { # loop from 2, as the 1st col is the outcome/response
    model = glm(y~ current_cohort[,k], family="binomial" )
    ps = getPvalues(model)
    pvals = c(pvals, ps[2])
  }
  pvals = c(1,pvals) # need to add a dummy p value for the intercept, just so that the indices from these univariate regression line up with the multireg case
  return(as.numeric(pvals))
}


##############################################################
# Loading & processing data
##############################################################

# Load misc data
geneNameLookup = read.table( geneNameLookupLoc  ,header=T)  # Load in the missenseSNP -> Gene Name lookup table
uc_vs_cd = read.table( uc_vs_cdLoc  ,header=F) # Load in the Data to separate UC vs CD cases (1=UC, 2=CD, anyone not on this list is a control)
annotationData = read.table( annotationDataLoc  ,header=F) # annotation file that stores information on the missense/eQTL SNPs (their IDs, effect sizes and effect alleles)

# depending on the mode in which we run, create an 'exclusion' list of people, we don't want in the analysis
if (subPhenoMode != 0 ) {
exclusionList = uc_vs_cd[ which(uc_vs_cd$V3 !=subPhenoMode) ,2] # this will create indi IDs to exlude either UC or CD cases
} else { exclusionList = c()}

# if we are running in 'subphenoMode' 3 (we are interested in running UC vs CD), create a list of people who have UC to be used as controls
uc_as_controls= c()
if (subPhenoMode == 3 ) {
  uc_as_controls = uc_vs_cd[ which(uc_vs_cd$V3 ==1) ,2]
}

#load Covariates file that which has the PCs and the batch effects
covsFile_orig =  read.table( covsFileLoc  ,header=T) 
if ( otherPCsAddedAsCovars  == TRUE) { covsFile = covsFile_orig } else {  covsFile = covsFile_orig[,1:7] } # exclude anything beyond the 7th col, which are the other PCs

# some hardcoded vars
cohorts = c("GWAS1","GWAS2","GWAS3","UKBB", "IBDSeq") 
tests = c("missense_SNP", "eQTL_SNP", "missense_x_eQTL", "#BAD_HAPLO")
tests_Chris = c("missense_hap", "eQTL_hap", "missense_x_eQTL")
cols_name = c("RSid", "test", "GWAS1"  ,  "GWAS2"  ,   "GWAS3"  ,  "UKBB"   , "IBDSeq"   , "Combined"  ,   "Meta", "significant")
sampleSizes = c(length=length(cohorts)) # empty array that will store the number of Cases in each cohort to be used as a weight for the Meta-p calculation
useCovs= TRUE # if covariates should be used (this is hardcoded to be true, but can disable it for debugging purposes)

# to be able to exclude the cohort columns when we run cohort specific regressions:
chris_cohort_idx = c(5:8) # the columns that are GWAS1, GWAS2, GWAS3 and UKBB, (IE the 4 cohort batch effect covars)
cohort_idx = c(6:9) # the indices for Chris' design matrix (different indices for the same terms as above)

# Setup the Results tables that will include results from all cohorts and all variants
results_all = cbind(NULL)
results_all_Chris = cbind(NULL)



##############################################################
# Main part of the script
##############################################################

subfolders = list.dirs(path =baseDir, full.names = F, recursive = F) # get a list of subfolders
for( i in 1:length(subfolders) ) {
  print(subfolders[i])
  missenseID =  subfolders[i] 
  tableLoc = paste(baseDir , missenseID , "/" , "_results_" , missenseID , "_TABLE", sep="" )
  tableLoc_Chris = paste(baseDir , missenseID , "/" , "_resultsChris_" , missenseID , "_TABLE", sep="" )
  
  # get the gene's name linked to the missense variant:
  geneName = toString( as.character(geneNameLookup[which(geneNameLookup$missenseSNP == missenseID),2] ) )
  # find missense var's chrom from the annotation data (as the eQTL SNP IDs are always chrom:coords ,we can get the chrom by just splitting the string...)
  annotDat = annotationData [which(annotationData$V1 == missenseID),4]
  chrom = as.numeric( unlist( strsplit(as.character(annotDat),":") )[1] )
  missenseID_display = paste(missenseID,"\n(chrom",chrom," / ",geneName,")", sep ="") # produce the name to be displayed in the results table, that should have the missenseSNP id, the chrom and the gene name
 
  # load data
  table =  read.table( tableLoc  ,header=T) 
  table_Chris =  read.table( tableLoc_Chris  ,header=T) 
  
  ###############
  # Process data: match the design matrices above, to the covariates file, by merging them by the individuals' IDs

  # I) #BAD_HAPLO model
  # need to load the .fam file with the IDs that is specific for current variant ( as each variant had a different missing rate exclusion rate)
  famFile =  read.table( paste(baseDir , missenseID , "/missense_eqtl_final.fam" , sep="" )  ,header=F) 
  table_withIDs = cbind(as.character(famFile$V2),table) # append the Indi IDs to the data, as the same .fam file was used for phasing the sample order should match
  
  # if we are running UC vs CD, we exclude all controls
  if (subPhenoMode == 3 ) { exclusionList = famFile[ which(famFile$V6 ==1) ,2]}
  table_withIDs = table_withIDs[!table_withIDs$`as.character(famFile$V2)` %in% exclusionList,]
  table_withIDs[table_withIDs$`as.character(famFile$V2)` %in% uc_as_controls,2] = 0 # set the UC people as controls for the UC vs CD runs ( the 'uc_as_controls' is normally empty for all other scenarios)
  
  # merge data table with covariates file based on the Individual IDs (to add the covariates back into the model)
  merged = merge(table_withIDs,covsFile,by.x="as.character(famFile$V2)",by.y="IID") 
  merged_cleaned = merged[,c(-1,-7)] # only keep the actual columns used for regression, exclude the IDs that are in column 1 and 7
  
  # II) Chris' haplo model: as each individual has 2 rows  and merge works with unique IDs, we can't do it in one go...
  # split the data, so that we have 1 haplo / indi, and merge it 2x with the covars then paste it back together
  table_Chris_hap1 = table_Chris[seq(1, nrow(table_Chris), 2), ] # take the 1st haplo for each indi
  table_Chris_hap2 = table_Chris[seq(2, nrow(table_Chris), 2), ] # take the 2nd haplo for each indi
  table_Chris_hap1_withIDs = cbind(as.character(famFile$V2),table_Chris_hap1) # add to both of them the indi IDs from the fam file
  table_Chris_hap2_withIDs = cbind(as.character(famFile$V2),table_Chris_hap2) 
  # if we are running UC vs CD, we exclude all controls
  table_Chris_hap1_withIDs = table_Chris_hap1_withIDs[!table_Chris_hap1_withIDs$`as.character(famFile$V2)` %in% exclusionList,]
  table_Chris_hap2_withIDs = table_Chris_hap2_withIDs[!table_Chris_hap2_withIDs$`as.character(famFile$V2)` %in% exclusionList,]
  table_Chris_hap1_withIDs[table_Chris_hap1_withIDs$`as.character(famFile$V2)` %in% uc_as_controls,2] = 0 # set the UC people as controls for the UC vs CD runs
  table_Chris_hap2_withIDs[table_Chris_hap2_withIDs$`as.character(famFile$V2)` %in% uc_as_controls,2] = 0 # set the UC people as controls for the UC vs CD runs
  
  # Merge the 2 haplos with the covariates separately
  merged_chris_hap1 = merge(table_Chris_hap1_withIDs,covsFile,by.x="as.character(famFile$V2)",by.y="IID") 
  merged_chris_hap2 = merge(table_Chris_hap2_withIDs,covsFile,by.x="as.character(famFile$V2)",by.y="IID") 
 
  # put the two back together, in a way so that they are interlaced IE hap2 follows hap1 for an individual
  merged_chris_hap1_2 = merged_chris_hap1[FALSE,] # create blank table in the right shape
  odd_rows = seq(1, nrow(table_Chris_hap1_withIDs)*2, 2)  # generate the odd and even rows as selection masks
  even_rows = seq(2, nrow(table_Chris_hap2_withIDs) *2, 2)
  merged_chris_hap1_2[ odd_rows,] = merged_chris_hap1[c(1:nrow(merged_chris_hap1)),] # paste in the 1st haplotype in the odd rows (with covariates)
  merged_chris_hap1_2[ even_rows,] = merged_chris_hap2[c(1:nrow(merged_chris_hap2)),] # paste in the 2nd haplotype in the even rows (with covariates)
  merged_Chris_cleaned =  merged_chris_hap1_2[,c(-1,-6)] # only keep the actual columns used for regression, exclude the IDs that are in column 1 and 6

  ###############################################
  
  ###############
  # Fit models for each cohort:  
  results_cohort <- data.frame(matrix(ncol = length(cohorts)+5, nrow = 4))   # +5 as =1 for RS id, +1 for the tests, and +1 for the 'overall' significant, +1 for Combined, and + 1 for meta
  results_cohort[,1] = rep(missenseID_display,4)
  results_cohort[,2] = tests
  results_cohort_Chris <- data.frame(matrix(ncol = length(cohorts)+5, nrow = 3))
  results_cohort_Chris[,1] = rep(missenseID_display,3)
  results_cohort_Chris[,2] = tests_Chris

  for( j in 1:length(cohorts) ) { # loop for each cohort, split the design matrices based on the cohort batch terms, IE to extract only people who are GWAS1, we get people who have a '1' in that column, and exclude everyone else
    print(cohorts[j]) 
    if (cohorts[j] == "IBDSeq") { # IBDseq does not have a term for it in the design matrix, so to get them by getting people that are '0' for all cohorts
      print("IBDSeq does not have a cohort, to avoid colinearity") 
      current_cohort = merged_cleaned[ which(merged_cleaned[ cohorts[1]]==0 & merged_cleaned[ cohorts[2]]==0 & merged_cleaned[ cohorts[3]]==0 & merged_cleaned[ cohorts[4]]==0),-cohort_idx] # -cohort_idx, as we don't want the cohort terms in the design matrix when we are only looking at a single cohort
      current_cohort_Chris = merged_Chris_cleaned[ which(merged_Chris_cleaned[ cohorts[1]]==0 & merged_Chris_cleaned[ cohorts[2]]==0 & merged_Chris_cleaned[ cohorts[3]]==0 & merged_Chris_cleaned[ cohorts[4]]==0),-chris_cohort_idx]
    } else {
    current_cohort = merged_cleaned[ which(merged_cleaned[cohorts[j]]==1),-cohort_idx]
    current_cohort_Chris = merged_Chris_cleaned[ which(merged_Chris_cleaned[cohorts[j]]==1),-chris_cohort_idx]
    }
   
    # Fit #BAD_HAPLO model for current cohort
    if (MARGINAL) { pvals = getMarginalPvals(current_cohort) } else { # if we are fitting marginal models
    myModel_glm_current = glm(y ~.,current_cohort, family="binomial" ) # if we fit the normal multireg GLM
    pvals = getPvalues(myModel_glm_current) # extract p-values from the model (need a separate function for it, as if some of the terms are NaN then the idiot R will return less p-values that potentially refer to the wrong terms...)
    }
    
    # same to fit Chris' haplo model for current cohort
    if (MARGINAL) { pvals_Chris = getMarginalPvals(current_cohort_Chris) } else {
    myModel_glm_current_chris = glm(y ~.,current_cohort_Chris, family="binomial" )
    ### ADD VARIANCE CORRECTION TO THIS MODEL HERE ###
    pvals_Chris = getPvalues(myModel_glm_current_chris)
    }
    
    # if there aren't enough cases (say 50), we replace the values with NAs (as sometimes we get p vals ~1 which would just mess up the meta analysis) (the 50 is used as GWAS1/2 which are 'pure' UC/CD cohorts still have a few of the other subtypes in them)
    if ( length ( which(current_cohort$y==1) ) < 50 || length ( which(current_cohort$y==0) ) < 50 ) {
      print( paste("not enough cases or controls for cohort:", cohorts[j]) )
      pvals = rep(NA, length(pvals))
      pvals_Chris = rep(NA, length(pvals_Chris))
    }
    
    # add the p-values for into the results from this cohort
    results_cohort[,(j+2)] = pvals[2:5]
    results_cohort_Chris[,(j+2)] = pvals_Chris[2:4]
    
    # if we are doing this the first time, the save number of cases for each cohort, as that will be used to weight the P-values
    if ( i == 1) { sampleSizes[j]= length ( which(current_cohort$y==1) )  }
  }
  
  # if we are not supposed to use Covariates for the combined test, then use the original design matrices 
  if ( useCovs == FALSE) {
    merged_cleaned = table
    merged_Chris_cleaned = table_Chris
  } 
  
  # Fit #BAD_HAPLO model for combined cohort
  if (MARGINAL) { pvals = getMarginalPvals(merged_cleaned) } else {
  myModel_glm = glm(y ~.,merged_cleaned, family="binomial" )
  pvals = getPvalues(myModel_glm)
  }  

  # same to fit Chris' haplo model for combined cohort
  if (MARGINAL) { pvals_Chris = getMarginalPvals(merged_Chris_cleaned) } else {
  myModel_glm_Chris = glm(y ~.,merged_Chris_cleaned, family="binomial" )
  ### ADD VARIANCE CORRECTION TO THIS MODEL HERE ###
  pvals_Chris = getPvalues(myModel_glm_Chris)
  }
  
  # add the combined results into the results
  results_cohort[,8] = pvals[2:5]
  results_cohort_Chris[,8] = pvals_Chris[2:4] 
  
  # get the meta p-values using Stouffer's) method # sum (w * z(p)) / sqrt(sum (w * w))
  for( j in 1:nrow(results_cohort) ) {
    pvals = as.numeric( results_cohort[j,3:7] )
    if ( sum(!is.na(pvals) ) < 2 ) { meta_p = NA } else {
      meta_p = as.numeric( sumz(pvals,sampleSizes, na.action =na.omit)[2] )
    }
    
    results_cohort[j,9] = meta_p
  }
  # same thing for Chris' model
  for( j in 1:nrow(results_cohort_Chris) ) {
    pvals = as.numeric( results_cohort_Chris[j,3:7] )
    if ( sum(!is.na(pvals) ) < 2 ) { meta_p = NA } else {
    meta_p = as.numeric( sumz(pvals,sampleSizes, na.action =na.omit)[2] )
    }
    results_cohort_Chris[j,9] = meta_p
  }
  
  # Cosmetics: if either the meta or the combined p-vals for the interaction terms were significant, ass a 1 into the last column (this will be used to format the table)
  if ( is.na(results_cohort[3,8]) == FALSE && results_cohort[3,8] < sigThreshold || is.na(results_cohort[4,8]) == FALSE && results_cohort[4,8] < sigThreshold || is.na(results_cohort[3,9]) == FALSE && results_cohort[3,9] < sigThreshold || is.na(results_cohort[4,9]) == FALSE && results_cohort[4,9]< sigThreshold )  {
    results_cohort[,10] = 1
  } else { results_cohort[,10] = 0 }

  if ( is.na(results_cohort_Chris[3,8]) == FALSE && results_cohort_Chris[3,8] < sigThreshold  || is.na(results_cohort_Chris[3,9]) == FALSE && results_cohort_Chris[3,9] < sigThreshold  )  {
    results_cohort_Chris[,10] = 1
  } else { results_cohort_Chris[,10] = 0 }
  

  # correct the estimates for Chris' model for variance inflation:
 
  results_all = rbind(results_all,results_cohort)
  results_all_Chris = rbind(results_all_Chris,results_cohort_Chris)
}


# write out the raw data
colnames(results_all) =cols_name
colnames(results_all_Chris) =cols_name
write.table(results_all, paste(outLoc,"resultsTable",sep=""), sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE) 
write.table(results_all_Chris, paste(outLoc,"resultsTable_Chris",sep=""), sep=" ", row.names = FALSE, col.names = FALSE, quote = FALSE) 
print(paste("written results to:", outLoc))


##############################################################
# Utility functions & Cosmetics
##############################################################

generateModelPredictions_BADHAPLO = function (current_cohort) { # generate model predictions for all genotypes 
  myModel_glm_current = glm(y ~.,current_cohort, family="binomial" ) # fit model
  datausedForPrediction = current_cohort
  modelForPrediction = myModel_glm_current
  summary(modelForPrediction)
  mean_PC1 = mean(datausedForPrediction$PC1)
  testIndiColnames = colnames(datausedForPrediction[1,])
  indinames = c()
  testindis = NULL
  
  for( i in 1:3 ) { # go through the missense genotypes (MM, Mm, ee)
    missenseGenotype = i -1 # 0,1 and 2
    for( j in 1:3 ) { # go through the eQTL genotypes (EE, Ee, ee)
      eQTLGenotype = j -1 # 0,1 and 2
      testIndi = c(0,missenseGenotype, eQTLGenotype, missenseGenotype * eQTLGenotype, 0 , mean_PC1 )
      testIndi_name = paste("M", missenseGenotype," / E",eQTLGenotype, sep="")
      indinames = c(indinames,testIndi_name )
      testindis = rbind(testindis, testIndi)
    }
  }
  
  rownames(testindis) = indinames
  testindis = as.data.frame((testindis))
  colnames(testindis) =  testIndiColnames
  predicted_y_all = predict(modelForPrediction, newdata = testindis, type = "response", se.fit = TRUE)
  predicyed_y_table = cbind(predicted_y_all$fit, predicted_y_all$se.fit)
  colnames(predicyed_y_table) = c("y_hat", "SE")

  write.table(predicyed_y_table, paste(outLoc,"predictions.csv", sep=""), sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE) 
  return(predicyed_y_table)
}


convertToScientific = function (res) { # converts numbers into scientific notation, and also casts them as character otherwise flextable ignores formatting
  results_scientific = res
  for( j in 3:(ncol(results_scientific) -1) ) {
    charCol = results_scientific[,j]
    charCol = formatC(as.numeric(charCol), format = "e", digits = 2)
    charCol = gsub("NA", "", charCol)
    results_scientific[,j] =  charCol
  }
  return(results_scientific)
}


displayTable = function (res, title, excludedCols = -1) { # creates a pretty HTML table using the flextable library
  colNames_used =  names(res)
  if  (excludedCols != -1 ) {
  colNames_used = colNames_used[-excludedCols]
  }
  
  hapFit <- flextable(res , col_keys =colNames_used) # , col_keys = cols_name[1:10]
  hapFit <- add_header_lines(hapFit,  values = title ) 
  hapFit

  hapFit <- theme_vanilla(hapFit)
  hapFit
  hapFit <- merge_v(hapFit, j = c("RSid" ) )
  hapFit <- autofit(hapFit)
  hapFit

  hapFit <- color(hapFit, i= ~ significant == 1, j = ~  RSid, color = "red")
  hapFit <-  bold(hapFit,i= ~ significant == 1,j = ~  RSid, bold = TRUE)

  hapFit <- color(hapFit, i= ~ as.numeric(Combined) < sigThreshold, j = ~  Combined, color = "red")
  hapFit <-  bold(hapFit,i= ~ as.numeric(Combined) < sigThreshold,j = ~  Combined, bold = TRUE)

  hapFit <- color(hapFit, i= ~ as.numeric(Meta) < sigThreshold, j = ~  Meta, color = "red")
  hapFit <-  bold(hapFit,i= ~  as.numeric(Meta) < sigThreshold,j = ~  Meta, bold = TRUE)
  hapFit <- align(hapFit, i =1, align = "center", part = "header")

  hapFit <- void(hapFit, j = "significant", part = "all" )

  hapFit
}


# build predictions for the main #BAD_HAPLO model
predicyed_y_table = generateModelPredictions_BADHAPLO(current_cohort)

# display results in a pretty table
displayTable(convertToScientific(results_all), "UC: #BAD_HAPLO model", -1)

##############################################################
# DEBUG: data for local testing
##############################################################

# covsFileLoc = "../data/misc/covarOut_covsAll"
# sigThreshold = 0.005
# uc_vs_cdLoc = "../data/misc/gwas1_3_ibdseq_ukbb_ukbb_merged_all_subtypes.phe"
# subPhenoMode = 1
# geneNameLookupLoc ="../data/misc/missense_gene_lookup.csv"
# MARGINAL=FALSE # if marginal/univariate regression should be fit instead of multiple regression
# otherPCsAddedAsCovars = FALSE # if the other PCs are to be used (other than PC1)
# 
# # the rs2476601 example for Chris
# baseDir = "../data/SNPs/"
# outLoc = "../results/"
# annotationDataLoc = "../data/misc/missense_eqtl_annotation.txt"
