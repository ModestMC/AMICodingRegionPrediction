## Usage
1. Configure input parameters as desired. These appear at the top of **predictCodingCrossSpecies.m**:
	- **seqLength** - Sequence length used for both training and testing the SVMs (default: 1000 basepairs)
	- **setSize** - The number of sequences drawn from the coding and noncoding parent sequences for both the training and test sets (default: 2000)
	- **minK/maxK** - Minimium and maximum lag *k* used to generated the AMI-derived profiles. The profiles consist of values for all lags between and including these extrema (default: 1/16).
	- **forceDownload** - If true, assembly files will be downloaded from NCBI even if the directory **NCBIGenomes** and the file **SpeciesList.csv** already exist, and the list in **TaxonomicIDList.csv** contains no IDs that do not occur in **SpeciesList.csv**. If false, assembly files will only be downloaded if one of those conditions is not met (default: false).
	- **taxIdFile** - File name of the taxonomic ID list (default: **TaxonomicIDList.csv**).
	- **outputFile** - File name of the output results (default: **CodingRegion_CrossSpecies.xlsx**).
2. Compile a list of taxonomic IDs corresponding to species of interest and place them in the taxonomic ID list file. The file consists of a single column, with one ID per row.
3. Run **predictCodingCrossSpecies.m** from the MATLAB command line
4. Examine results in output file **CodingRegion_CrossSpecies.xlsx**. Each sheet corresponds to a different profile (eaAMI, eAMI, or AMI) and classifier (SVM or Euclidean distance) pair. Each row provides performance metrics for a classifier trained on sequences drawn from the species indicated in the row header. Each set of 3 columns shows the results (AUC, sensitivity, and specificity) when testing that classifier on sequences drawn from each test species. If a species provided in the taxonomic ID list does not appear in the results file, it is likely that an assembly could not be found for it (indicated by the **validData** column of **SpeciesList.csv**).     

*Note*: If the configuration parameters are left at their defaults and **TaxonomicIDList.csv** is *not* modified, the script will use the assembly files provided in the **NCBIGenomes** directory. This is a collection of a few small genomes for species included in **TaxonomicIDList.csv**.  

## Example
Suppose one is interested in determining the effectiveness of cross-species coding region predictions for *S. cerevisiae S288C*, *S. enterica subsp. enterica serovar Typhimurium str. 14028S*, and *V. cholerae O1 str. C6706*. First, replace the list of IDs in **TaxonomicIDList.csv** with those for the species of interest (559292, 588858, and 1124478). Next, delete the **NCBIGenomes** directory and the file **SpeciesList.csv** (or set **forceDownload** to true). Lastly, run **predictCodingCrossSpecies.m**. If interested in the effects of the input parameters (e.g. sequence length), modify those parameters and re-run the script. If **forceDownload** is set to false, the script will skip the data compilation steps and proceed to training the classifiers. Note that the output file name should be changed to avoid overwriting previous results.  
