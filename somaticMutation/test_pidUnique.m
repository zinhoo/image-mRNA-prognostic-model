clear

load mutData.mat

pid = mutData.pid;

upid = unique(pid);

tcode = cellfun(@(x) x(14:15), pid, 'UniformOutput', false);
unique(tcode);