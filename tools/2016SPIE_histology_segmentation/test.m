clear
close all
tic

% I = imread('sample.jpg');
% I = imread('E:\KIRP_subimages_tissue/TCGA-2Z-A9JG-01A-01-TS1.08510C04-9A3C-451D-88A2-90F899A36A99_3.jpg');
% im = I(1:800, 1000:2000, :);

% I = imread('D:\matlab\kirp\survivalAnalysis\morphological\images_3\35fea_meanStdNuRegionProps_meanMinorAxisLength_high/TCGA-DZ-6133-01A-01-TS1.e666f310-f987-45dd-afa6-ed1e01f1e06b_1.jpg');
% im = I(1:800, 1000:2000, :);

% I = imread('E:\KIRP_subimages_tissue/TCGA-2Z-A9J6-01A-01-TS1.C50B97C0-C145-4939-BECD-587F9B482E01_0.jpg');
% I = imread('E:\KIRP_subimages_tissue/TCGA-5P-A9K3-01A-01-TS1.C3894CA5-A741-4BE4-9F65-175C36FB8A20_0.jpg');
I = imread('E:\KIRP_subimages_tissue/TCGA-AL-3467-01A-01-TS1.a4e25066-7012-4cb8-911b-3c317adb081f_1.jpg');

% imshow(im);

bw = hmt(I, [], true);

hold on

% stats = regionprops(bw, 'area', 'centroid');
% c = struct2cell(stats)';
% areas = cell2mat(c(:, 1));
% centroids = cell2mat(c(:, 2));
% 
% text(centroids(:, 1), centroids(:, 2), num2str(areas), 'color', 'w');

figure, imshow(bw);
figure, imshow(I);

% I = imread('sample.jpg');
% segmentedImage = hmt(I, [], true);
toc