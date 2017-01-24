% Sort the log-rank test results based on p value and add feature names.
clear

res = load('im_logrankRes.txt');
imFeasName = textread('../../imageFeatures/imFeasName.txt', '%s');

[v, ind] = sort(res(:, 2));

feaInd = res(ind, 1);
feaName = imFeasName(feaInd);
p = v;

im_logrankRes = table(feaInd, feaName, p);
save im_logrankRes im_logrankRes
