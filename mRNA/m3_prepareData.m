% Prepare data, remove non-solid tumor samples
clear

eiGene = load('Eigene_Normalized_SCC30_wnormalized.txt')';
strc = load('filtered_geneExp_KIRC.mat');
sid = strc.sampleID;

sid = cellfun(@(x) x(1:15), sid, 'uniformOutput', false);
codeTumor = cellfun(@(x) x(14:15), sid, 'uniformOutput', false);

% only retain solid tumor ('01')
ind = strcmp(codeTumor, '01');
pid = sid(ind);
eiGene = eiGene(ind, :);

rnaEigen = table(pid, eiGene);
save rnaEigen rnaEigen
