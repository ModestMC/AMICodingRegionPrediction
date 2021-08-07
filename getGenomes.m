function getGenomes(taxaID)
%GETGENOMES For each taxonomic ID in the input array, retrieve the genomic 
% FASTA and GTF file from NCBI

	folder = fullfile(fileparts(which(mfilename)), 'NCBIGenomes');

	species = table;
	species.taxaID = taxaID;
	species.speciesName  = cell(size(species.taxaID));
	species.speciesName  = cell(size(species.taxaID));
	species.organismName = cell(size(species.taxaID));
	species.kingdom      = cell(size(species.taxaID));
	species.accession    = cell(size(species.taxaID));
	species.validData    = false(size(species.taxaID));

	for i=1:length(species.taxaID)
		try
			fprintf('Obtaining genome for Taxonomic ID %d\n', species.taxaID(i));
			baseURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
			eutil = 'esearch.fcgi?';
			dbParam = 'db=genome';
			termParam = ['&term=txid', num2str(species.taxaID(i)), '[orgn]'];
			usehistoryParam = '&usehistory=y';
			esearchURL = [baseURL, eutil, dbParam, termParam, usehistoryParam];
			searchReport = webreadRetry(esearchURL);
			ncbi = regexp(searchReport,...
						  '<QueryKey>(?<QueryKey>\w+)</QueryKey>\s*<WebEnv>(?<WebEnv>\S+)</WebEnv>',...
						  'names');

			summary = webreadRetry([baseURL...
							   'esummary.fcgi?db=genome&rettype=gb&retmode=text&WebEnv=', ncbi.WebEnv,...
							   '&query_key=',ncbi.QueryKey]);

			genome = regexp(summary,...
							['<Item Name="Organism_Name" Type="String">(?<speciesName>.*?)</Item>.*',...
							 '<Item Name="Organism_Kingdom" Type="String">(?<kingdom>.*?)</Item>.*',...
							 '<Item Name="Assembly_Accession" Type="String">(?<accession>.*?)</Item>'],...
							'names');

			species.speciesName{i} = genome.speciesName;
			species.kingdom{i} = genome.kingdom;

			dbParam = 'db=assembly';
			filterStr = ' AND "latest refseq"[filter]';
			sortStr = '&sort=significance';
			esearchURL = [baseURL, eutil, dbParam, termParam, filterStr, sortStr, usehistoryParam];

			searchReport = webreadRetry(esearchURL);
			ncbi = regexp(searchReport,...
						  '<QueryKey>(?<QueryKey>\w+)</QueryKey>\s*<WebEnv>(?<WebEnv>\S+)</WebEnv>',...
						  'names');

			summary = webreadRetry([baseURL...
							   'esummary.fcgi?db=assembly&rettype=gb&retmode=text&WebEnv=', ncbi.WebEnv,...
							   '&query_key=',ncbi.QueryKey]);

			results = regexp(summary,...
							 ['<AssemblyAccession>(?<accession>.*?)</AssemblyAccession>.*?',...
							 '<Organism>(?<organism>.*?)</Organism>.*?',...
							 '<FtpPath_RefSeq>(?<ftpPath>.*?)</FtpPath_RefSeq>'],...
							 'names');

			result = results(1); 
			species.organismName{i} = result.organism;
			species.accession{i} = result.accession;
			
			result.ftpPath = strrep(result.ftpPath, 'ftp://ftp.ncbi.nlm.nih.gov/', '');
			
			ftpObj = ftp('ftp.ncbi.nlm.nih.gov');
			temp = split(result.ftpPath, '/');
			genomeFilename = [temp{end}, '_genomic.fna.gz'];
			gtfFilename = [temp{end}, '_genomic.gtf.gz'];

			mget(ftpObj, [result.ftpPath, '/', genomeFilename], folder);
			mget(ftpObj, [result.ftpPath, '/', gtfFilename], folder);
			
			copyfile(fullfile(folder, result.ftpPath, genomeFilename), folder);
			copyfile(fullfile(folder, result.ftpPath, gtfFilename), folder);
			
			species.validData(i) = true;
		catch
			warning('Data not found for taxonomic ID %d', species.taxaID(i));
		end
	end

	gunzip([folder, '/*.gz']);
	delete([folder, '/*.gz']);
	rmdir(fullfile(folder, 'genomes'), 's');

	writetable(species, 'SpeciesList.csv');
end


function text = webreadRetry(url)
    retries = 0;
    success = false;
    pause(0.1);
    
    while(~success && retries < 10)
        try
            text = webread(url);
            if(~contains(lower(text), 'error'))
                success = true;
            else
                error(text)
            end
        catch
            pause(0.5);
            retries = retries + 1;
            disp('Retrying...')
        end     
    end
    
    if (~success)
        error('Read failed');
    end
end
       
          
