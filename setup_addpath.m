% Add the tools to matlab search path

[a, b, c] = fileparts(mfilename('fullpath'));
root = a;

addpath(genpath(fullfile(root, 'tools', 'openslide-win64-20160612')));
addpath(genpath(fullfile(root, 'tools', 'fordanic-openslide-matlab-976ac90')));
addpath(genpath(fullfile(root, 'tools', '2016SPIE_histology_segmentation')));
addpath(fullfile(root, 'tools', 'heatmap'));

