% correlate survival-related image features with eigengenes
clear
close all

imFeasName = textread('../../imageFeatures/imFeasName.txt', '%s');

load('../../imRnaSM.mat')
imFeas = imRnaSM.imFeas;
eiGene = imRnaSM.eiGene;

logRes = load('../screenImFeatures/im_logrankRes.txt');
imFeasSur = imFeas(:, logRes(:, 1));
feaNamesSur = imFeasName(logRes(:, 1));

% create metagene names
nGene = size(eiGene, 2);
geneNames = cell(nGene, 1);
for i = 1:nGene
    geneNames{i} = ['eigengene', num2str(i)];
end
    

rho = corr(eiGene, imFeasSur, 'type', 'spearman');
% show heatmap
heatmap(rho, feaNamesSur, geneNames,  [],...
    'Colorbar', true, 'TickAngle', 45, 'ShowAllTicks', true);
% show heatmap with rho > 0.3
rho2 = rho;
rho2(abs(rho2)<=0.3) = 0;
figure
heatmap(rho2, feaNamesSur, geneNames,  [],...
    'Colorbar', true, 'TickAngle', 45, 'ShowAllTicks', true);

save spearman_rho rho