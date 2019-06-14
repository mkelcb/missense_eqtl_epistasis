# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 20:30:14 2019

@author: mk23
"""
#: OBJECTIVE: Parses haplotype files to produce design matrices for missense-eqtl interaction tests
# Input: 
# 1) a directory that stores the SHAPEIT3 phased haplo files in subdirectories
# 2) .fam file that stores case/control status 
# 3) an annotation file that stores information on the missense/eQTL SNPs (their IDs, effect sizes and effect alleles)
# Output:
# 1) two Design matrices for regression: one for the #BAD_HAPLO model and one for Chris' Haplo model

# Example command:
# python scripts/produceDesignMatrix-Chris.py --annotationFileLoc data/misc/missense_eqtl_annotation.txt --baseDir data/SNPs/

######################################################

import argparse
from collections import OrderedDict
import numpy as np
import os

def loadPhenos(famFileLoc) : # loads in a fam File and stores for each indi if they are case or control
    caseControlStatus = []
    indiData = []
    with open(famFileLoc, "r") as file:
        for i in file:
            famArray = i.rstrip().split()
            caseControlStatus.append( int( famArray[5] ) -1 ) # the 6th column in the fam is the case control status. 1 is control, 2 is case, but we want 0 for control and 1 for case (so later it can be used as an index in a 0 based array... )
            indiData.append( [famArray[0],famArray[1]] )
           
    return(caseControlStatus, indiData) # signature is array of len n, 0 if control, 1 if case


def loadAnnotationData(annotationFileLoc) : # parses the datafile that contains all annotations, will produce a dictionary {missenseID: data} that stores for each missensevar (key theSNPID), their effect allele [0], effect size[1], and a list of eQTLs[2]: SNPID[0] with effect Allele[1] and Zscore[2]
    annotationData = OrderedDict()          # Annotation file signature: missesnseSNPID[0], missense_effectAllele[1], missense_effectSize[2], eQTLSNPId[3], eQTL_effectAllele[4], eQTL_effectSize[5]
    with open(annotationFileLoc, "r") as file:
        for i in file:
            annotationArray = i.rstrip().split()
            if annotationArray[0] not in annotationData : # if we have not encountered this type of gene before, we create it
                missenseSNP = []
                missenseSNP.append(annotationArray[1]) # missense_effectAllele
                missenseSNP.append( float(annotationArray[2]) ) # missense_effectAllele
                missenseSNP.append( [] ) # list of eQTLs
                annotationData[ annotationArray[0] ] = missenseSNP
            else :
                missenseSNP = annotationData[ annotationArray[0] ]
                
            # once either fetched or created the missense SNP, we add to its eQTL list the current eQTL SNP
            missenseSNP[2].append( annotationArray[3] ) # eQTLSNPId
            missenseSNP[2].append( annotationArray[4] ) # eQTL_effectAllele
            missenseSNP[2].append( float(annotationArray[5]) ) # eQTL_effectSize
          
    return(annotationData)


# .haps file signature:
    # https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample
# 17 rsID A B 1 1   # first 5 are describing the data, and then each pair are the haplotypes for an individual
# 17 rsID 0 1
# "1 0" means that the first haplotype carries the B/second allele
def produceDesignMatrix(missense_hapsFileLoc,eqtl_hapsFileLoc, annotationData,famFileLoc) : # takes in a location to the haps file, the annotation annotationData, and the individuals' case/control status
    caseControlStatus, indiData = loadPhenos(famFileLoc)

    missenseHapsData = []
    with open(missense_hapsFileLoc, "r") as file:
        for i in file:
            SNP_array = i.rstrip().split()
            SNP_id = SNP_array[0]
            SNP_id = SNP_id.split("_")[0] #  17:57885142_G_A ->  17:57885142  ( as we dont have the alleles in the annotation data)
            SNP_array[0] = SNP_array[1] = SNP_id
            missenseHapsData = SNP_array
      
    missenseID = missenseHapsData[0]
    missenseData  = annotationData[missenseID]
    
    eQTL_SNP_array = []
    with open(eqtl_hapsFileLoc, "r") as file:
        for i in file:
            SNP_array = i.rstrip().split()
            SNP_id = SNP_array[0]
            SNP_id = SNP_id.split("_")[0] #  17:57885142_G_A ->  17:57885142  ( as we dont have the alleles in the annotation data)
            SNP_array[0] = SNP_array[1] = SNP_id
            eQTL_SNP_array = SNP_array

    # setup data with the correct shapes for results
    regressionData = np.zeros ( (len(caseControlStatus), 5 ) ) # [0], pheno, [1]: missense, [2] eqtl, [3]missense * eqtl, [4]: badhaplocount
    regressionData_Chris = np.zeros ( (len(caseControlStatus) * 2, 4 ) ) # [0], pheno, 2 outcomes for each indi, and 3 cols: missense, eQTL and their interaction
    indis_BAD_HAPLO = np.zeros( len(caseControlStatus), dtype = np.int ) # create empty array with length of the number of indis, that will store the number of bad haplotypes an indi has
    missense_genotypeData = np.zeros(len(caseControlStatus),dtype = np.int)
    missense_eQTL_interaction = np.zeros(len(caseControlStatus),dtype = np.int)
    eQTL_genotypeData = np.zeros(len(caseControlStatus),dtype = np.int)
    
    # set up the the variables
    missense_effectAllele = missenseData[0] # 'G'
    missense_effectEffectSize = float(missenseData[1]) # this is the OR, eg 0.884
    missense_allele_0 = missenseHapsData[3] # 'A'
    missense_allele_1 = missenseHapsData[4] # 'G'
    eQTL_SNPId = missenseData[2][0] # 1:161480649'
    eQTL_effectAllele = missenseData[2][1] # G
    eQTL_effectEffectSize = missenseData[2][2] # -6.6025 # this is the Zscore, its sign tells if 
    eQTL_allele_0 = eQTL_SNP_array[3] # G
    eQTL_allele_1 = eQTL_SNP_array[4] # A
    
    # need to decide what is the 'bad' missense allele, -> one that is related to increased OR,  (IE the alt allele if that increases OR, or the ref allele if the alt allele decreases OR (protective))
    if missense_allele_0 == missense_effectAllele and missense_effectEffectSize > 1 or missense_allele_1 == missense_effectAllele and missense_effectEffectSize < 1 :
        missenseALllele_BAD = 0 # this needs to be numeric, and correspond to the haplotype data, which are coded 0 and 1
        missense_allele_bad_nucleotide = missense_allele_0 # the string representation eg 'G'
    else  : 
        missenseALllele_BAD = 1
        missense_allele_bad_nucleotide = missense_allele_1
        
    # decide on the Odds ratio of the allele that increases risk
    if missense_allele_bad_nucleotide == missense_effectAllele : OR = missense_effectEffectSize
    else : OR = 1/missense_effectEffectSize # if it is the ref allele we flip the OR

    # Decide on which allele is up or down regulating the epxression of the linked missense variant, this depends on if the zscore is negative or positive
    if eQTL_allele_0 == eQTL_effectAllele and eQTL_effectEffectSize > 0 or  eQTL_allele_1 == eQTL_effectAllele and eQTL_effectEffectSize < 0:
        eqtlAllele_upReg = 0 # if the effect allele is allele_1, AND it has a positive sign, or if it is allele_2, but that has a negative sign, in this case allele_1 is associated with upregulating of the expressions
        eqtlAllele_upReg_nucleotide = eQTL_allele_0 # the string representation eg 'G'
    else : 
        eqtlAllele_upReg = 1 
        eqtlAllele_upReg_nucleotide = eQTL_allele_1

    # loop through haplotypes and produce the deisgn matrix
    indiCounter = 0
    for k in range(5, len(eQTL_SNP_array) -1 , 2) : # loop through the remainder after the first 5 elements, where each  pair is a new indi in the .fam file
        # get the 2 haplotypes for the missense
        missense_hap_1 = int(missenseHapsData[k]) # this is 0 or 1, need to map this to good or bad
        missense_hap_2 = int(missenseHapsData[k+1]) 
        missense_genotypeData[indiCounter] = missense_hap_1 + missense_hap_2
        # get the 2 haplotypes for the eQTL
        eQTL_hap_1 = int(eQTL_SNP_array[k])
        eQTL_hap_2 = int(eQTL_SNP_array[k+1])
        eQTL_genotypeData[indiCounter] = eQTL_hap_1 + eQTL_hap_2
        missense_eQTL_interaction[indiCounter] = missense_genotypeData[indiCounter] * eQTL_genotypeData[indiCounter]
        
        # BAD_HAPLO is defined as having the bad missense allele on the same chrom of an eQTL that is upregulating that gene
        if missense_hap_1 == missenseALllele_BAD and eQTL_hap_1 == eqtlAllele_upReg : indis_BAD_HAPLO[indiCounter] += 1  # this is a counter as we can have 0, 1 or 2 bad haplotypes
        if missense_hap_2 == missenseALllele_BAD and eQTL_hap_2 == eqtlAllele_upReg : indis_BAD_HAPLO[indiCounter] += 1  

        # populate the design matrix for Chris' model: we don't keep the original haplo encoding of 0=ref allele, 1=alt allele, but recode it as 0=good allele, 1=bad/upreg allele, although this probably doesn't make a difference in this case
        haploIndex = k -5 # 0 based index
        regressionData_Chris[haploIndex,0] = regressionData_Chris[haploIndex +1,0] = caseControlStatus[indiCounter] # we have the same phenotype for both haplos
        regressionData_Chris[haploIndex,1] = int(missense_hap_1 == missenseALllele_BAD) # we code 1 for the bad missense allele
        regressionData_Chris[haploIndex+1,1] = int(missense_hap_2 == missenseALllele_BAD) # second haploytype for missense
        regressionData_Chris[haploIndex,2] = int(eQTL_hap_1 == eqtlAllele_upReg)  # code 1 for the upregulating eQTL allele
        regressionData_Chris[haploIndex+1,2] = int(eQTL_hap_2 == eqtlAllele_upReg) # secon haploytype for the eqTL
        # the interaction term, multiplying the two columns missenes * eQTL
        regressionData_Chris[haploIndex,3] = regressionData_Chris[haploIndex,1] * regressionData_Chris[haploIndex,2]
        regressionData_Chris[haploIndex+1,3] = regressionData_Chris[haploIndex+1,1] * regressionData_Chris[haploIndex+1,2]   
        indiCounter += 1
 
    # populate the design matrix for the genotype + #BAD_HAPLO model
    regressionData[:,0] =caseControlStatus
    regressionData[:,1] =missense_genotypeData
    regressionData[:,2] =eQTL_genotypeData
    regressionData[:,3] =missense_eQTL_interaction
    regressionData[:,4] =indis_BAD_HAPLO
    
    # diagnostics vars
    badHaploCountTotal = np.sum(indis_BAD_HAPLO)
    missense_2nd_alleleFreq = np.sum(missense_genotypeData) /  (len(caseControlStatus) *2 ) # the haplos basically count the 2nd allele's frequency (which is NOT necessarily the major allele), it is usually actually the minor allele
    eqtl_2nd_alleleFreq = np.sum(eQTL_genotypeData) /  (len(caseControlStatus) *2 )
    
    # decide on the allele freqs of the effect/upreg alleles
    missense_effect_AF = missense_2nd_alleleFreq
    if missenseALllele_BAD == 0 : missense_effect_AF = 1-missense_2nd_alleleFreq # if the bad missense allele is 0, then missense_2nd_alleleFreq, is the good allele's freq, so to get the bad one, it will be 1-good_freq
    
    eqtl_upreg_AF = eqtl_2nd_alleleFreq
    if eqtlAllele_upReg == 0 : eqtl_upreg_AF = 1-eqtl_2nd_alleleFreq
    
    # produce some diagnostiscs for a sanity check
    diagnostics = "Allele freq of missense var ("+missenseID+") effect allele (" + missense_allele_bad_nucleotide+ ")" + " with OR: " + str( np.round(OR,3) ) + " is: " + str( np.round(missense_effect_AF,2) ) +"\n"
    diagnostics += "Allele freq of eqtl var ("+eQTL_SNPId+") upregulating allele (" + eqtlAllele_upReg_nucleotide+ ") is: " + str( np.round(eqtl_upreg_AF,2) ) +"\n"
    bad_HaploFreq =  badHaploCountTotal / (len(caseControlStatus) *2 )
    diagnostics += "Bad haplotype arrangement frequency: " + str( np.round( bad_HaploFreq , 2) )
    print(diagnostics)
    diagnostics_data = [ [missenseID, missense_allele_bad_nucleotide,  missense_effect_AF ] , [eQTL_SNPId, eqtlAllele_upReg_nucleotide, eqtl_upreg_AF] , ["BAD Haplo feq exp/obs", missense_effect_AF * eqtl_upreg_AF, bad_HaploFreq]]
    
    return(regressionData, regressionData_Chris, diagnostics, diagnostics_data)

     
# writes out the design matrix to be used for the genotype + #BAD_HAPLO regression, as well as some diagnostics
def resultsToDisk(resultsFileLoc, missenseResults, diagnostics, diagnostics_data) :
    with open(resultsFileLoc + diagnostics_data[0][0] +"_TABLE", "w") as file: 
        file.write("y" + "\t" + "missenseSNP" + "\t" + "eQTLSNP" + "\t" + "missense_x_eQTL" + "\t" + "BAD_HAPLO_count" + "\n")
        for j in range(len(missenseResults)) :
            file.write(str(missenseResults[j][0]) + "\t" + str(missenseResults[j][1]) + "\t" + str(missenseResults[j][2]) + "\t" + str(missenseResults[j][3]) + "\t" + str(missenseResults[j][4]) + "\n")
       
    print("written out regression table for, haplo", diagnostics_data[0][0], " to:", resultsFileLoc + diagnostics_data[0][0] +"_TABLE")      

    # write out diagnostics too
    with open(resultsFileLoc + diagnostics_data[0][0] +"_diagnostics", "w") as file: 
            file.write(diagnostics)
            
    with open(resultsFileLoc + diagnostics_data[0][0] +"_diagtable", "w") as file: 
        for j in range(len(diagnostics_data)) :
            file.write(str(diagnostics_data[j][0]) + "\t" + str(diagnostics_data[j][1]) + "\t" + str(diagnostics_data[j][2]) + "\n")


# writes out the design matrix to be used for Chris' Haplotype based 
def resultsToDisk_Chris(resultsFileLoc, missenseResults_Chris, missenseSNPName) :
    with open(resultsFileLoc + missenseSNPName +"_TABLE", "w") as file: 
        file.write("y" + "\t" + "missense_SNP" + "\t" + "eQTL_SNP" + "\t" + "missense_x_eQTL"  + "\n")
        for j in range(len(missenseResults_Chris)) :
            file.write(str(missenseResults_Chris[j][0]) + "\t" + str(missenseResults_Chris[j][1]) + "\t" + str(missenseResults_Chris[j][2]) + "\t" + str(missenseResults_Chris[j][3])  + "\n")
       
    print("written out regression table for Chris' model  to:", resultsFileLoc + missenseSNPName +"_TABLE")      


def runMain(args) :
    print('Producing Design matrices from haplotype data started:')  
    baseDir = args.baseDir
    annotationFileLoc = args.annotationFileLoc
    
    annotationData = loadAnnotationData(annotationFileLoc)
    for x in os.listdir(baseDir): # loop through all subfolders in the base directory, and produce design matrices for each
        print(x)
        missenseResults, missenseResults_Chris, diagnostics, diagnostics_data = produceDesignMatrix(baseDir +x +"/_missenseExtract.haps" ,baseDir +x +"/_eqtlExtract.haps" , annotationData, baseDir +x +"/"+"missense_eqtl_final.fam") 
        resultsToDisk(baseDir +x +"/"+"_results_", missenseResults, diagnostics, diagnostics_data)
        resultsToDisk_Chris(baseDir +x +"/"+"_resultsChris_", missenseResults_Chris, diagnostics_data[0][0])

if __name__ == '__main__':   
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotationFileLoc",required=True, help='The location of the annotation file that stores for each missense/eQTL SNP pair their effect size and effect alleles')  
    parser.add_argument("--baseDir", help='The directory that should contain a subfolder')  
    parser.set_defaults(func=runMain)
        
    args = parser.parse_args()
    args.func(args)



############################################################################
############################################################################

# DEBUG data to be loaded in locally
#famFileLoc = 'C:/softwares/Cluster/META_GWAS/gwas_chr17clean.fam'
   
## eQTLGen
#annotationFileLoc = 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/0eQTL-Missense_Epistasis/data/annotation/missense_eqtl_annotation.txt'
#baseDir = "C:/softwares/Cluster/missense_eqtl_epistasis/data/" 
  
# BLUEPRINT Monoctyes
#baseDir = "C:/softwares/Cluster/missense_eqtl_epistasis/data_bliueprint_mono/"  
#annotationFileLoc = 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/0eQTL-Missense_Epistasis/data/annotation/blueprint_mono/missense_eqtl_annotation.txt'

## BLUEPRINT neutro
#baseDir = "C:/softwares/Cluster/missense_eqtl_epistasis/data_bliueprint_neutro/"  
#annotationFileLoc = 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/0eQTL-Missense_Epistasis/data/annotation/blueprint_neutro/missense_eqtl_annotation.txt'

## BLUEPRINT T cells
#baseDir = "C:/softwares/Cluster/missense_eqtl_epistasis/data_bliueprint_tcell/"  
#annotationFileLoc = 'C:/Users/mk23/GoogleDrive_phd/PHD/Project/0eQTL-Missense_Epistasis/data/annotation/blueprint_tcell/missense_eqtl_annotation.txt'

# picking out an individual variant
#x = "rs6025"
#missense_hapsFileLoc = baseDir +x +"/_missenseExtract.haps"
#eqtl_hapsFileLoc= baseDir +x +"/_eqtlExtract.haps"

