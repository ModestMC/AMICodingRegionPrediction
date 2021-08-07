%PREDICTCODING For the taxonomic ID specified:
% 1. [If Necessary] Download genome from NCBI
% 2. [If Necessary] Compile coding and noncoding parent sequences
% 3. Generate AMI profiles for a random selection of training sequences 
%    drawn from the coding and noncoding parent sequences
% 4. Train an SVM on the profiles in the training set
% 5. Use that SVM to predict whether each sequence in the specified FASTA
%    file is coding
% 6. Output prediction scores to specified output file (CSV, XLSX, or TXT) 


%% Input Parameters
seqLength = 1000;
setSize = 2000;
minK = 1;
maxK = 16;
forceDownload = false;
trainTaxaID = 400667; % Acinetobacter baumannii
inputFile = 'GCF_000046845.1_ASM4684v1_cds_from_genomic.fna'; % Acinetobacter baylyi
outputFile = 'CodingRegion_PredictionScores.csv';


%% Get Data (if necessary)
genomeFolder = 'NCBIGenomes';
speciesListFile = 'SpeciesList.csv';

if (isfile(speciesListFile))
    species = readtable(speciesListFile);
    if (~ismember(trainTaxaID, species.taxaID))
        forceDownload = true;
    end
end

if (forceDownload || ~isfolder(genomeFolder) || ~isfile(speciesListFile))
	getGenomes(trainTaxaID);
    compileCodingNoncoding;
end

species = readtable(speciesListFile);

tIdx = find(species.taxaID == trainTaxaID);
if (isempty(tIdx) || ~species.validData(tIdx))
    error('Could not find assembly for species with taxonomic ID %d', trainTaxaID);
end

if (~isfile(fullfile(genomeFolder, [strrep(species.accession{tIdx}, 'A', 'F'), '_coding_noncoding.fasta'])))
	getGenomes(trainTaxaID);
    compileCodingNoncoding;
    
    if (~isfile(fullfile(genomeFolder, [strrep(species.accession{tIdx}, 'A', 'F'), '_coding_noncoding.fasta'])))
        error('Could not find coding/noncoding file for species with taxonomic ID %d', trainTaxaID);
    end
end

if (~isfile(inputFile))
    error('Could not find input file "%s"', inputFile);
end


%% Setup
profiles = {'eaAMI', 'eAMI', 'AMI'};

profileValsTrain = cell(length(profiles));
profileValsTest = cell(length(profiles));
cSVM = cell(length(profiles));


%% Train SVMs
disp(['Training SVMs for ' species.speciesName{tIdx}]);
[headers, sequences] = fastaread(fullfile(genomeFolder, [strrep(species.accession{tIdx}, 'A', 'F'), '_coding_noncoding.fasta']));

allClass = [-1*ones(setSize, 1); ones(setSize, 1)];

coding = sequences{strcmp(headers, 'Coding')};
noncoding = sequences{strcmp(headers, 'Noncoding')};

trainSequence = cell(length(allClass)/2, 2);

for i=1:length(trainSequence)
    idx = randi(length(noncoding) - seqLength(l));
    trainSequence{i, 1} = noncoding(idx:idx+seqLength(l)-1);
    idx = randi(length(coding) - seqLength(l));
    trainSequence{i, 2} = coding(idx:idx+seqLength(l)-1);
end

trainSequence = reshape(trainSequence, 1, [])';

groupSequences = true;
if (seqLength(l) >= 800)
    groupSequences = false;
end

for t=1:length(profiles)
    [profileValsTrain{t}, AMIList] = getAMI(trainSequence, minK, maxK, profiles{t}, groupSequences);

    svmArgList = {'KernelFunction', 'linear', 'standardize', true, 'KernelScale', 'auto', 'BoxConstraint', .1 };
    cSVM{t} = fitcsvm(profileValsTrain{t}, allClass, svmArgList{:});
end


%% Predict Input Sequences
[headers, sequences] = fastaread(inputFile);
results = table(cellfun(@length, sequences)', 'VariableNames', {'Sequence Length (bp)'}, 'RowNames', headers');
for t=1:length(profiles)
    disp(['Making predictions using SVM on ', profiles{t}, ' for ' species.speciesName{tIdx}]);
    
    [profileVals, ~] = getAMI(sequences, minK, maxK, profiles{t}, false);
    [~, scores] = predict(cSVM{t}, profileVals);
    results.([profiles{t}, ' SVM Scores']) = scores(:,2);
end


%% Export results
if (endsWith(outputFile, '.xlsx', 'IgnoreCase', true))
    options = {'WriteMode', 'replacefile'};
elseif (endsWith(outputFile, '.csv', 'IgnoreCase', true))
    options = {'WriteMode', 'overwrite'};
else
    options = {'WriteMode', 'overwrite', 'Delimiter', '\t'};
end
writetable(results, outputFile, 'WriteRowNames', true, 'WriteVariableNames', true, options{:});
