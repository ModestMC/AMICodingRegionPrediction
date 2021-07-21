function compileCodingNoncoding()
%COMPILECODINGNONCODING Compile coding and noncoding parent sequences for
% each species listed in 'SpeciesList.csv'

    genomeFolder = 'NCBIGenomes';
    speciesListFile = 'SpeciesList.csv';

    species = readtable(fullfile(fileparts(which(mfilename)), 'SpeciesList.csv'));
    files = {dir(genomeFolder).name};

    for i=1:length(species.accession)
        if (species.validData(i))
            try
                fprintf('Compiling sequences for %s\n', species.speciesName{i});
                accessionGCF = strrep(species.accession{i}, 'A', 'F');
                genomicFile = files(contains(files, accessionGCF) & endsWith(files, 'genomic.fna'));
                gtfFile = files(contains(files, accessionGCF) & contains(files, '.gtf'));

                gtfAnnot = GFFAnnotation(fullfile(genomeFolder, gtfFile{1}));
                gtfStruct = getData(gtfAnnot);
                gtfStruct = gtfStruct(startsWith({gtfStruct.Feature}, 'CDS'));

                [headers, sequences] = fastaread(fullfile(genomeFolder, genomicFile{1}));
                genome = upper(sequences);

                if (~iscell(sequences))
                    headers = {headers};
                    genome = upper({sequences});
                end

                nonMit = ~contains(headers, 'mitochond', 'IgnoreCase', true);
                headers = headers(nonMit);
                genome = genome(nonMit);

                coding = {};
                noncoding = {};

                c = 1;
                found = 0;
                chromLimit = length(headers);
                if (dir(fullfile(genomeFolder, genomicFile{1})).bytes > 100e6)
                    chromLimit = min(length(headers), 3);
                end

                cdsAnnot = [];
                while (c <= length(headers) && found < chromLimit)
                    cAcc = extractBefore(headers{c}, ' ');
                    gtfStructChrom = gtfStruct(startsWith({gtfStruct.Reference}, cAcc));
                    [~, idx] = sort(cell2mat({gtfStructChrom.Start}));
                    gtfStructChrom = gtfStructChrom(idx);

                    if (~isempty(gtfStructChrom))
                        found = found + 1;
                        isCoding = false(length(genome{c}), 1);

                        for j=1:length(gtfStructChrom)
                            isCoding(gtfStructChrom(j).Start:gtfStructChrom(j).Stop) = true;
                            temp = genome{c}(gtfStructChrom(j).Start:gtfStructChrom(j).Stop);
                            if (gtfStructChrom(j).Strand == '-')
                                temp = seqrcomplement(temp);
                            end

                            coding = [coding, temp];
                        end

                        noncoding = [noncoding, genome{c}(~isCoding)];
                    end
                    c = c + 1;
                end

                noncoding =strjoin(noncoding, '');
                coding = strjoin(coding, '');
                noncoding = [noncoding, seqrcomplement(noncoding)];

                fastaName = [accessionGCF, '_coding_noncoding.fasta'];
                if (isfile(fullfile(genomeFolder, fastaName)))
                    delete(fullfile(genomeFolder, fastaName));
                end
                fastawrite(fullfile(genomeFolder, fastaName), {'Coding', 'Noncoding'}, {coding, noncoding});

            catch
               warning('Could not find %s', species.speciesName{i});
            end
        end
    end    
end
