# AMICodingRegionPrediction
Predict coding regions using SVMs on AMI-derived profiles of test sequences
 
## About
Running the script **predictCodingCrossSpecies.m** does the following for each taxonomic ID in **TaxonomicIDList.csv**:
1. [*If Necessary*] Download genome assemblies from NCBI and write summary information to **SpeciesList.csv** (by calling **getGenomes.m**)
2. [*If Necessary*] Compile coding and noncoding parent sequences (by calling **compileCodingNoncoding.m**)
3. Generate AMI profiles for a random selection of training sequences  drawn from the coding and noncoding parent sequences (by calling **getAMI**)
4. Train an SVM on the profiles in the training set
5. Use that SVM to predict whether test sequences drawn from each of the other species is coding or noncoding
6. Output prediction metrics to **CodingRegion_CrossSpecies.xlsx**

Running the script **predictCoding.m** does the following for the taxonomic ID specified:
1. [If Necessary] Download genome from NCBI (by calling **getGenomes.m**)
2. [If Necessary] Compile coding and noncoding parent sequences (by calling **compileCodingNoncoding.m**)
3. Generate AMI profiles for a random selection of training sequences drawn from the coding and noncoding parent sequences (by calling **getAMI**)
4. Train an SVM on the profiles in the training set
5. Use that SVM to predict whether each sequence in the specified multi-FASTA file is coding
6. Output prediction scores to specified output file (CSV, XLSX, or TXT) 

## Requirements
* MATLAB R2020a or later
* Internet connection (if downloading files from NCBI)

## License
MIT
