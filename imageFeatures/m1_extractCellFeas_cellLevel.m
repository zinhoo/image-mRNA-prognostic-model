% Extract 10 cell features for each image in 'imageInfo'
clear
tic

% nucleus area in real size
mpp_40 = 0.2525; % 40X
realMin = 80*mpp_40^2;
realMax = 1100*mpp_40^2;
realPSize = 2000*mpp_40;


% Load openslide library
openslide_load_library();

strc = load('../imageAndClinicalInfo/imageInfo.mat');
imInfo = strc.imageInfo;
dirData = '../imageAndClinicalInfo/images/';

% process each image
nFile = size(imInfo, 1);
for i = 1:nFile
    t1 = tic;
    file = imInfo.file{i};
    parts = regexp(file, '/', 'split');
    filename = parts{2}(1:end-4);
    
    % test if the image is already processed
    if exist(['extractCellFeas_cellLevel/', filename, '.mat'], 'file');        
        fprintf('image %d/%d already processed\n', i, nFile);
        continue;
    end
    
    slidePtr = openslide_open([dirData, file]);
    mpp = imInfo.mppX(i);
    width = imInfo.width(i);
    height = imInfo.height(i);
    ps = round(realPSize/mpp);
    
    %% select patches using thumnail
    % Get thumbnail
    thumnail = openslide_read_associated_image(slidePtr, 'thumbnail');
    thumnail = thumnail(:, :, 2:4);
    [ht, wt, ~] = size(thumnail);
    ratio = width/wt;
    
    % patch size in thumbnail
    pst = round(ps/ratio);

    % select patches of interest in thumbnail image
    [X, Y] = meshgrid(1:pst:wt-pst+1, 1:pst:ht-pst+1);
    xy = [X(:), Y(:)];
    d1 = size(xy, 1);
    indOk = zeros(d1, 1);
    for j = 1:d1
        r = xy(j, 2);
        c = xy(j, 1);
        rows = r:r+pst-1;
        cols = c:c+pst-1;
        tile = thumnail(rows, cols, :);
        tMean = mean(tile, 3);
        if sum(tMean(:)>210) < pst*pst/2
            indOk(j) = 1;
        end
    end
    xy = xy(indOk==1, :)-1; % move upper-left point to (0, 0)
    xy = round(xy*ratio);
    
    % remvoe coordinates exceeding width or height
    indOk = xy(:, 1)<=width-ps & xy(:, 2)<=height-ps;
    xy = xy(indOk, :);
    
    %% extract cell features for each tile
    d1xy = size(xy, 1);
    cellFeas = cell(d1xy, 1);
    for ixy = 1:d1xy
        x = xy(ixy, 1);
        y = xy(ixy, 2);
        tile = openslide_read_region(slidePtr, x, y, ps, ps, 0);
        tile = tile(:, :, 2:4);
        areaMin = round(realMin/mpp^2);
        bw = hmt(tile, areaMin);

        statsR = regionprops('table', bw, tile(:, :, 1), 'area', 'MeanIntensity',...
            'MajorAxisLength', 'MinorAxisLength', 'centroid');
        statsR.Area = statsR.Area*mpp^2;
        statsR.MajorAxisLength = statsR.MajorAxisLength*mpp;
        statsR.MinorAxisLength = statsR.MinorAxisLength*mpp;
        
        indOk = statsR.Area <= realMax;
        statsR = statsR(indOk, :);
        if size(statsR, 1) < 80
            continue;
        end
        
        statsG = regionprops('table', bw, tile(:, :, 2), 'MeanIntensity');
        statsB = regionprops('table', bw, tile(:, :, 3), 'MeanIntensity');
        
        statsG = statsG(indOk, :);
        statsB = statsB(indOk, :);

        feas1_7 = [statsR.Area, statsR.MajorAxisLength, statsR.MinorAxisLength,...
            statsR.MajorAxisLength ./ statsR.MinorAxisLength,...
            statsR.MeanIntensity, statsG.MeanIntensity, statsB.MeanIntensity];


        % last 3 features derived from Delaunay graph
        centroids = statsR.Centroid;
        feas8_10 = zeros(size(statsR, 1), 3);
        DT = delaunayTriangulation(double(centroids));
        E = edges(DT);
        for k = 1:size(centroids, 1)
            edgesCell = E(sum(E==k, 2)~=0, :);
            dist = zeros(size(edgesCell, 1), 1);
            for m = 1:numel(dist)
                p1 = centroids(edgesCell(m, 1), :);
                p2 = centroids(edgesCell(m, 2), :);
                dist(m) = norm(p1-p2);
            end
            feas8_10(k, :) = [mean(dist), max(dist), min(dist)];
        end

        cellFeas{ixy, 1} = [feas1_7, feas8_10*mpp];        
    end
        
    cellFeas(cellfun(@isempty, cellFeas)) = []; % remove empty elements
    save(['extractCellFeas_cellLevel', '/', filename, '.mat'], 'cellFeas');
    fprintf('image %d/%d processed, time %f\n', i, nFile, toc(t1));
end


toc

% Close whole-slide image, note that the slidePtr must be removed manually
openslide_close(slidePtr)
clear slidePtr

% Unload library
openslide_unload_library