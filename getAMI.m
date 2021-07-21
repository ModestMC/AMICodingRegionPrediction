function [AMIProfile, AMIList] = getAMI(sequences, minK, maxK, type, group)
%GETAMI Generate AMI profiles for a given set of sequences

    numPairs = 4^2;
    numElements = (maxK - minK + 1);
    
    seqLen = cellfun(@length, sequences);
    maxSeqLen = max(seqLen);
    
    if (strcmp(type, 'eaAMI') || strcmp(type, 'eAMI'))
    	numElements = numElements * numPairs;
    end
            
    
    if (sum(iscell(sequences)) == 0)
        sequences = { sequences };
    end
        
    AMIProfile = zeros(length(sequences), numElements);
    
    AMIList = cell(numElements, 1);
    
    
    seqNum = cell(length(sequences), 1);
    marginalProbs = zeros(length(sequences), 4);
    
	for i=1:length(sequences)
        seq = upper(sequences{i});
        seqNum{i} = zeros(1, length(seq));

        seqNum{i}(seq == 'T') = 1;
        seqNum{i}(seq == 'C') = 2;
        seqNum{i}(seq == 'G') = 3;
        
        for nt=1:size(marginalProbs, 2)
        	marginalProbs(i, nt) = sum(seqNum{i} == nt-1) / length(seqNum{i});  
        end
    end
    
    if (nargin > 4)
        if (group == true)
            temp = -1*4^(maxK+1)*ones(length(sequences), maxSeqLen);

            for i=1:length(sequences)
                temp(i, 1:seqLen(i)) = seqNum{i};   
            end

            seqNum = temp;
            clear temp;
            
            width = maxSeqLen - maxK;
            fullConvMatrix = 4 * eye(maxSeqLen);
            fullConvMatrix = fullConvMatrix(:, 1:width);
        
            fullConvMatrix(minK+1:maxSeqLen+1:end) = 1;
        end
    end
    

        
	for k=minK:maxK
        
        if (strcmp(type, 'AMI'))
            AMIList{k} = sprintf('%i', k);
        end
        
        for kmer=1:numPairs
            kmerNum = dec2base(kmer-1, 4, 2);

            kmerSeq = kmerNum;

            kmerSeq(kmerNum == '0') = 'A';
            kmerSeq(kmerNum == '1') = 'T';
            kmerSeq(kmerNum == '2') = 'C';
            kmerSeq(kmerNum == '3') = 'G';

            if (strcmp(type, 'eaAMI') || strcmp(type, 'eAMI'))
                AMIList{ numPairs * (k-minK) + kmer } = sprintf('k=%i-%s', k, kmerSeq);
            end
        end
        
        multiplierVect = zeros(k+1, 1);
        multiplierVect(1) = 1;
        multiplierVect(end) = 4;
        
        if ((nargin > 4 && group == false) || nargin <= 4)
            for i=1:length(sequences)
                AMILen = length(seqNum{i}) - maxK;
                
                AMIIndices = conv(multiplierVect, seqNum{i});
                AMIIndices = AMIIndices(k+1:k+AMILen);

                AMIThisK = zeros(1, numPairs);

                for kmer=1:numPairs
                    AMIThisK(kmer) = sum(AMIIndices == kmer-1);  
                end

                AMIThisK = AMIThisK ./ AMILen;


                if (strcmp(type, 'eaAMI') || strcmp(type, 'AMI'))
                    AMIThisK = AMIThisK .* log2(AMIThisK ./ kron(marginalProbs(i, :), ones(1, 4)) ./ repmat(marginalProbs(i, :), 1, 4));
                    AMIThisK(isnan(AMIThisK)) = 0;
                end

                if (strcmp(type, 'eaAMI') || strcmp(type, 'eAMI'))
                    AMIProfile(i, numPairs * (k-minK)+1:numPairs * (k-minK+1)) = AMIThisK;
                else
                    AMIProfile(i, k+1-minK) = sum(AMIThisK);
                end
            end
        else
            AMILen = seqLen - maxK;

        	fullConvMatrixa = convmtx(multiplierVect', maxSeqLen);
            fullConvMatrixa = fullConvMatrixa(:, k+1:maxSeqLen-maxK+k);

            AMIIndices = seqNum * fullConvMatrix;

            fullConvMatrix(k+1:maxSeqLen+1:end) = 0;
            fullConvMatrix(k+2:maxSeqLen+1:end) = 1;

            AMIThisK = zeros(length(sequences), numPairs);

            for kmer=1:numPairs
                AMIThisK(:, kmer) = sum(AMIIndices == kmer-1, 2);  
            end

            AMIThisK = AMIThisK ./ repmat(AMILen, 1, numPairs);


            if (strcmp(type, 'eaAMI') || strcmp(type, 'AMI'))
                AMIThisK = AMIThisK .* log2(AMIThisK ./ (kron(marginalProbs, ones(1, 4)) .* repmat(marginalProbs, 1, 4)));
                AMIThisK(isnan(AMIThisK)) = 0;
            end

            if (strcmp(type, 'eaAMI') || strcmp(type, 'eAMI'))
                AMIProfile(:, numPairs * (k-minK)+1:numPairs * (k-minK+1)) = AMIThisK;
            else
                AMIProfile(:, k+1-minK) = sum(AMIThisK, 2);
            end
        end
    end
end

