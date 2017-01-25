% Sort the log-rank test results based on p value
clear

res = load('logrankRes.txt');
[v, ind] = sort(res(:, 2));

feaInd = res(ind, 1);
p = v;

logrankRes = table(feaInd, p);
save logrankRes logrankRes
