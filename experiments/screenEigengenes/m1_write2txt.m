% Write eigengene data into .txt file which will be processed in R.
% For each feature, patients are divided into two groups by the median of
% the feature.
clear

load('../../imRnaSM.mat');
eiGene = imRnaSM.eiGene;
time = imRnaSM.cliInfo.time/30;
death = imRnaSM.cliInfo.death;

label = zeros(size(eiGene));
for i = 1:size(eiGene, 2)
    fea = eiGene(:, i);
    cutoff = prctile(fea, 50);
    label(fea<cutoff, i) = 1;
    label(fea>=cutoff, i) = 2;
end

dlmwrite('rdata.txt', [time, death, label], '\t');