% find high correlation pairs (correlation > 0.3)
clear

load spearman_rho
logrank = load('../screenImFeatures/im_logrankRes.txt');
imFeasName = textread('../../imageFeatures/imFeasName.txt', '%s');

[rs, cs] = find(abs(rho)>0.3);
ind = find(abs(rho)>0.3);
val = rho(ind);

[geneImSorted, idx] = sortrows([rs, cs]);

highCorr = [geneImSorted, val(idx)];

% map index to original image feature index
highCorr(:, 2) = logrank(highCorr(:, 2), 1);
name = imFeasName(highCorr(:, 2));

metaGene = highCorr(:, 1);
imFea = name;
corr = highCorr(:, 3);
highCorr_tab = table(metaGene, imFea, corr);
save highCorr_tab highCorr_tab 
% writetable(tab, 'highCorr.txt', 'delimiter', 'tab');