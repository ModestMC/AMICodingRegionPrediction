## Usage

#### predictCodingCrossSpecies.m
1. Configure input parameters as desired. These appear at the top of **predictCodingCrossSpecies.m**:
	- **seqLength** - Sequence length used for both training and testing the SVMs (default: 1000 basepairs)
	- **setSize** - The number of sequences drawn from the coding and noncoding parent sequences for both the training and test sets (default: 2000)
	- **minK/maxK** - Minimum and maximum lag *k* used to generated the AMI-derived profiles. The profiles consist of values for all lags between and including these extrema (default: 1/16).
	- **forceDownload** - If true, assembly files will be downloaded from NCBI even if the directory **NCBIGenomes** and the file **SpeciesList.csv** already exist, and the list in **TaxonomicIDList.csv** contains no IDs that do not occur in **SpeciesList.csv**. If false, assembly files will only be downloaded if one of those conditions is not met (default: false).
	- **taxIdFile** - File name of the taxonomic ID list (default: **TaxonomicIDList.csv**).
	- **outputFile** - File name of the output results. Must be XLSX (default: **CodingRegion_CrossSpecies.xlsx**).
2. Compile a list of taxonomic IDs corresponding to species of interest and place them in the taxonomic ID list file. The file consists of a single column, with one ID per row.
3. Run **predictCodingCrossSpecies.m** from the MATLAB command line
4. Examine results in output file. Each sheet corresponds to a different profile (eaAMI, eAMI, or AMI) and classifier (SVM or Euclidean distance) pair. Each row provides performance metrics for a classifier trained on sequences drawn from the species indicated in the row header. Each set of 3 columns shows the results (AUC, sensitivity, and specificity) when testing that classifier on sequences drawn from each test species. If a species provided in the taxonomic ID list does not appear in the results file, it is likely that an assembly could not be found for it (indicated by the **validData** column of **SpeciesList.csv**).     

*Note*: If the configuration parameters are left at their defaults and **TaxonomicIDList.csv** is *not* modified, the script will use the assembly files provided in the **NCBIGenomes** directory. This is a collection of a few small genomes for species included in **TaxonomicIDList.csv**.  

#### predictCoding.m
1. Configure input parameters as desired. These appear at the top of **predictCoding.m**. These are the same as for **predictCodingCrossSpecies.m** except for the following:
	- **trainTaxaID** - Taxonomic ID of species used for training the SVMs (default: 400667).
	- **inputFile** - Name of the multi-FASTA file containing sequences of interest to be scored using the SVMs (default: **GCF_000046845.1_ASM4684v1_cds_from_genomic.fna**).
	- **outputFile** - File name of the output results. Can be XLSX, CSV, or TXT (default: **CodingRegion_PredictionScores.xlsx**).
2. Run **predictCoding.m** from the MATLAB command line
3. Examine results in output file. The first column contains each of the headers from the multi-FASTA file. The second column contains the length of each sequence, in basepairs. Columns 3-5 contain the SVM scores for each sequence using the profile noted in the column header. Positive scores indicate that the sequence is predicted to be coding. Higher magnitude scores imply a higher probability that the sequence is coding or noncoding. Typically, eAMI will provide better results for short sequences (less than 250 basepairs), while eaAMI will provide better results for long sequences. 
      

## Examples

#### predictCodingCrossSpecies.m
Suppose one is interested in determining the effectiveness of cross-species coding region predictions for *S. cerevisiae S288C*, *S. enterica subsp. enterica serovar Typhimurium str. 14028S*, and *V. cholerae O1 str. C6706*. Suppose also we want the length of the sequences to be 200 bp and we want to use 1000 sequences for training and testing.  The input parameter setup is the first portion of the MATLAB program **predictCodingCrossSpecies.m**.  On line 13 we set **seqLength** to 200, and on line 14 we set **setSize** to 1000.

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
	Obtaining genome for Taxonomic ID 559292
	Obtaining genome for Taxonomic ID 588858
	Obtaining genome for Taxonomic ID 1124478
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

If we now wish to see what happens when we change some of the parameters, such as the sequence length, we can do that and rerun the program.  Because the sequences have already been downloaded and the corresponding information compiled, the program goes directly to training the SVM and generating the output file, which speeds up the process. Be aware that unless you change the name of the output file by setting **outputFile** to a new name the program will overwrite the previous output file.  If you want to force the program to re-download the files and recompile the associated information you can do that by either deleting the **NCBIGenomes** directory and the file **SpeciesList.csv** or by setting **forceDownload** to true.


#### predictCoding.m
Suppose one is interested in identifying potential genes in a newly-sequenced bacterial genome. A reasonable first step would be to scan the genome for long ORFs and compile these ORFs into a multi-FASTA file (e.g. **newGenomeORFs.fna**). We could then use **predictCoding.m** to predict whether each ORF is a protein-coding gene. We would need to determine a suitable species to be used for training. A reasonable choice would be the most closely related species for which an assembly and annotation exists in the NCBI database. Suppose that species is *Streptococcus pneumoniae D39* (taxonomic ID 373153).  

The input parameter setup is the first portion of the MATLAB program **predictCoding.m**. We would set **inputFile** to **newGenomeORFs.fna**, and **trainTaxaID** to 373153. If we want the output to be an XLSX file instead of CSV, we would set **outputFile** to, for example, **newGenomeORFScores.xlsx**.  A reasonable choice for **seqLength** is the median length of the sequences in the input multi-FASTA file. If this is unknown, one can run the program once, which will provide the length of each sequence in the output file. 2000 is likely a suitable choice for **setSize**, but a larger number may be better for a large training genome.        

Assuming this is the first time we are using this program we can simply run **predictCoding.m**.  The program will use the function **getGenomes.m** to download the training genome from NCBI.  In this example the command window in MATLAB will show the progression of the program:

	>> predictCoding
	Obtaining genome for Taxonomic ID 373153
	Compiling sequences for Streptococcus pneumoniae
	Training SVMs for Streptococcus pneumoniae
	Making predictions using SVM on eaAMI for Streptococcus pneumoniae
	Making predictions using SVM on eAMI for Streptococcus pneumoniae
	Making predictions using SVM on AMI for Streptococcus pneumoniae

The program will generate a spreadsheet (default name **CodingRegion_PredictionScores.csv**) with the prediction scores for the different profiles in each column. If we want to identify a few ORFs that are highly likely to be protein-coding genes, we could select the ORFs assigned the highest eaAMI SVM scores.  

If we now wish to see what happens when we change some of the parameters, such as the sequence length, we can do that and rerun the program. For example, if we are particularly interested in short ORFs, selecting a shorter training sequence length may produce more accurate results for those sequences. Because the sequences have already been downloaded and the corresponding information compiled, the program goes directly to training the SVMs and generating the output file. Be aware that unless you change the name of the output file by setting **outputFile** to a new name the program will overwrite the previous output file. If you want to force the program to re-download the files and recompile the associated information you can do that by either deleting the **NCBIGenomes** directory and the file **SpeciesList.csv** or by setting **forceDownload** to true.
