% Intersection of image, RNAseq, SM, and clinical data.

clear

load('../imageFeatures/imData.mat')
load('../mRNA/rnaEigen.mat')
load('../mutation/mutData.mat')
load('../imageAndClinicalInfo/cliInfo.mat')

% intersecct image, RNA, clinical data, and mutation data
[pid1, ind1, ind2] = intersect(imData.pid, rnaEigen.pid);
imFeas = imData.imFeas(ind1, :);
eiGene = rnaEigen.eiGene(ind2, :);

[pid2, ind1, ind2] = intersect(pid1, cliInfo.pid);
imFeas = imFeas(ind1, :);
eiGene = eiGene(ind1, :);
cliInfo = cliInfo(ind2, :);

[pid3, ind1, ind2] = intersect(pid2, mutData.pid);
imFeas = imFeas(ind1, :);
eiGene = eiGene(ind1, :);
cliInfo = cliInfo(ind1, :);
mut = mutData.mut(ind2, :);

imRnaSM.cliInfo = cliInfo;
imRnaSM.imFeas = imFeas;
imRnaSM.eiGene = eiGene;
imRnaSM.mut = mut;
imRnaSM.geneName = mutData.geneName;

save('imRnaSM', 'imRnaSM');
