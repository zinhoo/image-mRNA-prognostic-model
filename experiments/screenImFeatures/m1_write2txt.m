% Write images features into .txt file which will be processed in R
% For each feature, patients are divided into two groups by the median of
% the feature.
clear

load('../../imRnaSM.mat');
imFeas = imRnaSM.imFeas;
time = imRnaSM.cliInfo.time/30;
death = imRnaSM.cliInfo.death;

label = zeros(size(imFeas));
for i = 1:size(imFeas, 2)
    fea = imFeas(:, i);
    cutoff = prctile(fea, 50);
    label(fea<cutoff, i) = 1; % low group
    label(fea>=cutoff, i) = 2; % high group
end

dlmwrite('rdata.txt', [time, death, label], '\t');