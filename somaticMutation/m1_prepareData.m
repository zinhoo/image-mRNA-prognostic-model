% prepare mutation data
clear

tab = readtable('KIRC_coerced_oncotator_mutation_table.txt',...
    'ReadVariableNames', false, 'delimiter', 'tab');
[s1, s2] = size(tab);

pid = tab{1, 2:end}';
for i = 1:numel(pid)
    pid{i} = strrep(pid{i}, '_', '-');
end
    
geneName = tab{2:end, 1};
mut = zeros(s2-1, s1-1);

for i = 2:s2
    ind1 = strcmp(tab{2:end, i}, '1');
    mut(i-1, ind1) = 1;
    mut(i-1, ~ind1) = 0;
end

mutData.pid = pid;
mutData.geneName = geneName;
mutData.mut = mut;

save mutData mutData





