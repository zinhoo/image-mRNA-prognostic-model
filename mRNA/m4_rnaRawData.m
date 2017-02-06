% Prepare raw RNA data, remove non-solid tumor samples
tic
clear

strc = load('filtered_geneExp_KIRC.mat');
sid = strc.sampleID;
geneName = strc.GeneSymbol_f2;
geneExp = strc.geneExp_f2';

% remove '|' in geneName
for i = 1:numel(geneName)
    c = regexp(geneName{i}, '\|', 'split');
    geneName{i} = c{1};
end

pid = cellfun(@(x) x(1:15), sid, 'uniformOutput', false);
codeTumor = cellfun(@(x) x(14:15), sid, 'uniformOutput', false);

% only retain solid tumor ('01')
ind = strcmp(codeTumor, '01');
rnaRaw.pid = pid(ind);
rnaRaw.geneExp = geneExp(ind, :);
rnaRaw.geneName = geneName;

save rnaRaw rnaRaw
toc