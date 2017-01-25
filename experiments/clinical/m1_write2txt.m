% Write clinical data into .txt file which will be processed in R
clear

load('../../imRnaSM.mat');

time = imRnaSM.cliInfo.time/30;
death = imRnaSM.cliInfo.death;
grade = imRnaSM.cliInfo.grade;
stage = imRnaSM.cliInfo.stage;

tab = table(time, death, grade, stage);
writetable(tab, 'rdata.txt', 'WriteVariableNames', false);