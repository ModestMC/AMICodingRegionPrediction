# AMICodingRegionPrediction

Predict coding regions using SVMs on AMI-derived profiles of test sequences
 
## Usage
Running the script **predictCodingCrossSpecies.m** does the following for each taxonomic ID in **TaxonomicIDList.csv**:
1. [If Necessary] Download genomes from NCBI
2. [If Necessary] Compile coding and noncoding parent sequences
3. Generate AMI profiles for a random selection of training sequences  drawn from the coding and noncoding parent sequences
4. Train an SVM on the profiles in the training set
5. Use that SVM to predict whether test sequences drawn from each of the  other species is coding or noncoding
6. Output prediction metrics to **CodingRegion_CrossSpecies.xlsx**

Configurable input parameters appear at the top of the script. A limited number of genomic FASTA files are included. To download genomes for the remaining species in **TaxonomicIDList.csv** (or if additional taxonomic IDs are added), set **forceDownload** to true or delete the **NCBIGenomes** directory.

## License

MIT