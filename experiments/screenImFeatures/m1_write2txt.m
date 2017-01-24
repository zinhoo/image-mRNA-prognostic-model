% Write data into .txt file which will be processed in R
clear

load('../../imRnaSM.mat');
imFeas = imRnaSM.imFeas;
time = imRnaSM.cliInfo.time/30;
death = imRnaSM.cliInfo.death;

label = zeros(size(imFeas));
for i = 1:size(imFeas, 2)
    fea = imFeas(:, i);
    cutoff = prctile(fea, 50);
    label(fea<cutoff, i) = 1;
    label(fea>=cutoff, i) = 2;
end

dlmwrite('rdata.txt', [time, death, label], '\t');