% Prepare data for uni- and multivariate Cox analysis
clear

rnaName = {'CSNK2A1'; 'SPP1'; 'DEFB1'; 'PECAM1'; 'EDNRB'; 'TSPAN7'};
mutName = {'VHL'; 'PBRM1'; 'BAP1'; 'SETD2'; 'TP53'};

load('../../mRNA/rnaRaw.mat');
load('../../imRnaSM.mat');
rnaGeneName = rnaRaw.geneName;
mutGeneName = imRnaSM.geneName;

time = imRnaSM.cliInfo.time/30;
death = imRnaSM.cliInfo.death;

% rna raw
[pset, ia, ib] = intersect(rnaRaw.pid, imRnaSM.cliInfo.pid);
nr = numel(rnaName);
indRna = zeros(nr, 1);
for i = 1:nr
    indRna(i) = find(strcmp(rnaGeneName, rnaName{i}));
end
rnaSel = rnaRaw.geneExp(ia, indRna);
rnaSelBi = zeros(size(rnaSel));
% binarize
for i = 1:size(rnaSel, 2)
    med = median(rnaSel(:, i));
    rnaSelBi(rnaSel(:, i)<med, i) = 1;
    rnaSelBi(rnaSel(:, i)>=med, i) = 2;
end

% mut
nm = numel(mutName);
indMut = zeros(nm, 1);
for i = 1:nm
    indMut(i, 1) = find(strcmp(mutGeneName, mutName{i}));
end
mutSel = imRnaSM.mut(:, indMut);

% lasso-Cox grouping
coxGroup = load('../lassoCox/group_loo_imGene_related_norm.txt');

% grade12 vs. 34
grade = zeros(size(coxGroup));
ind34 = strcmp(imRnaSM.cliInfo.grade, 'G3') | strcmp(imRnaSM.cliInfo.grade, 'G4');
grade(ind34) = 1;

% stage12 vs. 34
stage = zeros(size(coxGroup));
ind34 = strcmp(imRnaSM.cliInfo.stage, 'Stage III') | strcmp(imRnaSM.cliInfo.stage, 'Stage IV');
stage(ind34) = 1;

dlmwrite('rdata.txt', [time, death, grade, stage, coxGroup, rnaSelBi, mutSel], '\t');

