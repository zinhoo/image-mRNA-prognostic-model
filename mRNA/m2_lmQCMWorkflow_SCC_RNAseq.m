% this code will take the curated matrix data after filtering out low varinance and low mean data,
% and search for high density clusters based on lmQCM algorithm, normalized first, then compute eigen
% genes and clusters, merge small clusters, and output merged cluster genes
% and eigen values. Remember to set gamma, t, lamda, beta and minimum
% clusters before run. 

%%%%%%%%%%% Step 1 - Compute PCC matrix %%%%%%%%%%%%%%%%

%%%% if there is NaN or missing data in the data matrix, use below %%%%%%
tic

% nGene=11611; %change here for gene list length
[sortData, sortInd] = sort(geneExp_f2'); %change here to your own matrix variable
[sortData, sampleData] = sort(sortInd);
sampleData = sampleData';

%%%% calculate the correlation matrix and setting the diagonal to be zeros
cMatrix = corr(sampleData', 'type', 'Spearman');
cMatrix(1 : size(cMatrix, 1) + 1 : end) = 0; 

%%%%%%%% if weigth normalization is needed, use below %%%%%%%%%%%%%%%%%%%%
cMatrix = abs(cMatrix);
D = sum(cMatrix);
D_half = 1./sqrt (D);
for i = 1 : size(cMatrix, 1)
    cMatrix(i, :) = cMatrix(i, :) * D_half(i);
end;
for i = 1 : size(cMatrix, 1)
    cMatrix(:, i) = cMatrix(:, i) * D_half(i);
end;

%%%%%%%% Step 2 - identify co-expression modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Algorithm parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 0.3; t = 1; lambda = 1;
%%%%%%%% Run the algorithm
C = localMaximumQCM(abs(cMatrix), gamma, t, lambda);

%%%%%%%% Step 3 - Merge the overlapped networks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Allowed overlap ratio threshold for merging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.4; minClusterSize = 10;

%%%%%%%% Sort and iteratively merge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizeC = zeros(1, length(C));
for i = 1 : length(C)
    sizeC(i) = length(C{i});
end;

[sortC, sortInd] = sort(sizeC, 'descend');
C = C(sortInd);

ind = find(sortC >= minClusterSize);

mergedCluster = C(ind);
mergeOccur = 1; 
currentInd = 0;

while mergeOccur == 1
    mergeOccur = 0;
    while currentInd < length(mergedCluster)
        currentInd = currentInd + 1;
        excludeInd = [];
        if (currentInd < length(mergedCluster))
            keepInd = 1 : currentInd;
            for j = currentInd+1 : length(mergedCluster)
                interCluster = intersect(mergedCluster{currentInd}, mergedCluster{j});
                if length(interCluster) >= beta*min(length(mergedCluster{j}), length(mergedCluster{currentInd}))
                    mergedCluster{currentInd} = union(mergedCluster{currentInd}, mergedCluster{j});
                    mergeOccur = 1;
                else
                    keepInd = [keepInd, j];
                end;
            end;
            mergedCluster = mergedCluster(keepInd);
            length(mergedCluster)
        end;
    end;
    sizeMergedCluster = zeros(1, length(mergedCluster));
    for i = 1 : length(mergedCluster)
        sizeMergedCluster(i) = length(mergedCluster{i});
    end;
    [sortSize, sortMergedInd] = sort(sizeMergedCluster, 'descend');
    mergedCluster = mergedCluster(sortMergedInd);
    currentInd = 0;
end;


%%%% Saving the gene modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputName = strcat('MergedGeneIDList_Normalized_SCC', num2str(round(gamma*100)), '_wnormalized.txt');
fid = fopen(outputName, 'w');
for i = 1 : length(mergedCluster)
    tmp = mergedCluster{i};
    for j = 1 : length(tmp)
        fprintf(fid, '%s\t', GeneSymbol_f2{tmp(j)}); %change this variable to your own gene list
    end;
    fprintf(fid, '\n');
end;
fclose(fid);

%%%% calculate and save the eigengenes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outputName = strcat('Eigene_Normalized_SCC', num2str(round(gamma*100)), '_wnormalized.txt');
fid = fopen(outputName, 'w');
for i = 1 : length(mergedCluster)
    tmp = mergedCluster{i};
    tmpData = geneExp_f2(tmp, :);  %change here to your own matrix variable name
    meanTmpData = mean(tmpData,1);
    centeredTmpData = tmpData - repmat(meanTmpData, size(tmpData,1), 1);
    for j = 1 : size(tmpData,1)
        centeredTmpData(j, :) = centeredTmpData(j, :)/norm(centeredTmpData(j, :));
    end;
    
    [u, s, v] = svd(centeredTmpData);
    signS = sign(corr(tmpData',v(:,1)));
    
    if (sum(signS) < 0)
        tmpEigenGene = -v(:,1);
    else
        tmpEigenGene = v(:,1);
    end

    for j = 1 : size(tmpData, 2)
        fprintf(fid, '%f\t', tmpEigenGene(j));
    end;
    fprintf(fid, '\n');
end;
fclose(fid);

%save the index for each cluster, maybe needed to fetch other values later
outputName = strcat('Cluster_indices_SCC', num2str(round(gamma*100)), '_wnormalized.txt');
fid = fopen(outputName, 'w');
for i = 1 : length(mergedCluster)
    tmp = mergedCluster{i};
    for j = 1 : length(tmp)
        fprintf(fid, '%d\t', tmp(j)); 
    end;
    fprintf(fid, '\n');
end;
fclose(fid);

toc

