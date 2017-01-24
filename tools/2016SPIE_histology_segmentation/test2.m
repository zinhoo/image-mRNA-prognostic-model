clear
close all
tic



I = imread('D:\matlab\kirp\survivalAnalysis\morFeas_v3\images_3\35fea_meanStdNuRegionProps_meanMinorAxisLength_high/TCGA-DZ-6133-01A-01-TS1.e666f310-f987-45dd-afa6-ed1e01f1e06b_1.jpg');
imLarge = I(1:800, 1000:2000, :);

I = imread('D:\matlab\kirp\survivalAnalysis\morFeas_v3\images_3\35fea_meanStdNuRegionProps_meanMinorAxisLength_low/TCGA-B3-3926-01A-01-TS1.970e79b4-ab56-49a8-820b-1a3e63c46056_2.jpg');
imSmall = I(1:800, 1000:2000, :);

bwLarge = hmt(imLarge, [], true);
hold on
stats = regionprops(bwLarge, 'area', 'centroid');
c = struct2cell(stats)';
areas = cell2mat(c(:, 1));
centroids = cell2mat(c(:, 2));
text(centroids(:, 1), centroids(:, 2), num2str(areas), 'color', 'w');

figure
bwSmall = hmt(imSmall, [], true);
hold on
stats = regionprops(bwSmall, 'area', 'centroid');
c = struct2cell(stats)';
areas = cell2mat(c(:, 1));
centroids = cell2mat(c(:, 2));
text(centroids(:, 1), centroids(:, 2), num2str(areas), 'color', 'w');

toc