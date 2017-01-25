% Write survival-associated features into .txt file which will be input to
% R for lasso-regularized Cox regression
% Features are standardized
clear

strc1 = load('../../imRnaSM.mat');
strc2 = load('../screenImFeatures/im_logrankRes.mat');
strc3 = load('../screenEigengenes/logrankRes.mat');

imFeas = strc1.imRnaSM.imFeas;
eiGene = strc1.imRnaSM.eiGene;

indIm = strc2.im_logrankRes.feaInd;
indGene = strc3.logrankRes.feaInd;

imFeas2 = zscore(imFeas(:, indIm));
eiGene2 = zscore(eiGene(:, indGene));
imGene = [imFeas2, eiGene2];

time = strc1.imRnaSM.cliInfo.time/30;
death = strc1.imRnaSM.cliInfo.death;

dlmwrite('rdata_imGene_related_norm.txt', [time, death, imGene], '\t');