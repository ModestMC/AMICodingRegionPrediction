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
Suppose one is interested in determining the effectiveness of cross-species coding region predictions for *S. cerevisiae S288C*, *S. enterica subsp. enterica serovar Typhimurium str. 14028S*, and *V. cholerae O1 str. C6706*. Suppose also we want the length of the sequences to be 200 bp and we want to use 1000 sequences for training and testing.  The input parameter setup is the first porion of the MATLAB program **predictCodingCrossSpecies.m**.  On line 13 we set **seqLength** to 200, and on line 14 we set **setSize** to 1000.

We now need the taxonomic ID for the sequences we are interested in.  We can get them from the NCBI taxonomy page.  If we do that we get the following set of taxonomic IDs

*S. cerevisiae S288C*: 559292  
*S. enterica subsp. enterica serovar Typhimurium str. 14028S*: 588858  
*V. cholerae O1 str. C6706*: 1124478  

We store this list of IDs in a .csv file.  The default name for the file is **TaxonomicIDList.csv**.  We can change the default name if we want to by changing the value of the variable **taxIdFile** on line 23.  In our case the file **TaxonomicIDList.csv** looks like this: 

taxaID  
559292  
588858  
1124478  

Assuming this is the first time we are using this program we can simply run **predictCodingCrossSpecies.m**.  The program will use the function **getGenomes.m** to download the genomes from NCBI.  In this example the command window in MATLAB will show the progression of the program:

	>> predictCodingCrossSpecies
	Obtaining genome for Taxonomoic ID 559292
	Obtaining genome for Taxonomoic ID 588858
	Obtaining genome for Taxonomoic ID 1124478
	Compiling sequences for Saccharomyces cerevisiae
	Compiling sequences for Salmonella enterica
	Compiling sequences for Vibrio cholerae
	Training SVM for Saccharomyces cerevisiae
	Training SVM for Salmonella enterica
	Training SVM for Vibrio cholerae
	Making predictions using SVM on eaAMI for Saccharomyces cerevisiae
	Making predictions using SVM on eaAMI for Salmonella enterica
	Making predictions using SVM on eaAMI for Vibrio cholerae
	Making predictions using SVM on eAMI for Saccharomyces cerevisiae
	Making predictions using SVM on eAMI for Salmonella enterica
	Making predictions using SVM on eAMI for Vibrio cholerae
	Making predictions using SVM on AMI for Saccharomyces cerevisiae
	Making predictions using SVM on AMI for Salmonella enterica
	Making predictions using SVM on AMI for Vibrio cholerae

The program will generate a spreadsheet (default name **CodingRegion_CrossSpecies**) with the results for the different profiles in different tabs.  It also generates a file called **SpeciesList.csv** with the species name, organism name, kingdom, and accession.

If we now wish to see what happens when we change some of the parameters such as the sequence length we can do that and rerun the program.  Because the sequences have already been downloaded and the corresponding information compiled the program goes directly to training the SVM and generating the output file speeding up the process. Be aware that unless you change the name of the output file by setting **outputFile** to a new name the program will overwrite the previous output file.  If you want to force the program to re-download the files and recompile the associated information you can do that by either deleting the **NCBIGenomes** directory and the file **SpeciesList.csv** or by setting **forceDownload** to true.
