% OPENSLIDE_EXAMPLE An example of how to use the MATLAB bindings in this library
% for working with whole-slide images

% Copyright (c) 2014 Daniel Forsberg
% danne.forsberg@outlook.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


clear
close all

% Get folder of this file
[folder,~,~] = fileparts(mfilename('fullpath'));

% Load openslide library
openslide_load_library();

% Set file name of whole-slide image to work with
% WSI = [folder,filesep,'test_data',filesep,'CMU-1-Small-Region.svs'];

% jpeg2000
WSI = 'E:\KIRP\nationwidechildrens.org_KIRP.diagnostic_images.Level_1.299.3.0\TCGA-A4-A5DU-01Z-00-DX1.ACCD40B6-CC28-4E89-B9A9-9199F0541075.svs';
% WSI = 'C:\Users\cheng\Desktop\type1\TCGA-2Z-A9J9-01Z-00-DX1.D6AA9708-04ED-477F-9FAF-B70F0DD8C1BE.svs';
% jpeg
% WSI = 'C:\Users\cheng\Desktop\type2\TCGA-Q2-A5QZ-01Z-00-DX1.66131B27-D89A-4BCC-A25A-1A4B9682199C.svs';

% Display vendor of whole-slide image
disp(['Vendor: ',openslide_detect_vendor(WSI)])

% Open whole-slide image
slidePtr = openslide_open(WSI);

% Get whole-slide image properties
[mppX, mppY, width, height, numberOfLevels, ...
    downsampleFactors, objectivePower] = openslide_get_slide_properties(slidePtr);

% Display properties
disp(['mppX: ',num2str(mppX)])
disp(['mppY: ',num2str(mppY)])
disp(['width: ',num2str(width)])
disp(['height: ',num2str(height)])
disp(['number of levels: ',num2str(numberOfLevels)])
% disp(['downsample factors: ',num2str(downsampleFactors)])
disp(['objective power: ',num2str(objectivePower)])

% Read a part of the image
ARGB0 = openslide_read_region(slidePtr,7308, 4917, 500, 500, 0);
ARGB1 = openslide_read_region(slidePtr,1780, 1029, 256, 256, 1);
ARGB2 = openslide_read_region(slidePtr,812, 512, 256, 256, 2);

white = openslide_read_region(slidePtr, 4409, 32506, 256, 256, 0);

imargb = openslide_read_region(slidePtr,7308, 4917, 224, 224, 0);
figure
imshow(imargb(:, :, 2:4));

% Display RGB part
figure
imshow(ARGB0(:,:,2:4))
figure
imshow(ARGB1(:,:,2:4))
figure
imshow(ARGB2(:,:,2:4))


% Get property names and display them
propertyNames = openslide_get_property_names(slidePtr);
disp(propertyNames(:))

% Get a property value
propertyValue = openslide_get_property_value(slidePtr,propertyNames{21});
disp(propertyValue)

propertyValue = openslide_get_property_value(slidePtr,propertyNames{40});
disp(propertyValue)

% Get associated images
associatedImages = openslide_get_associated_image_names(slidePtr);
disp(associatedImages(:))

% Get label image
label = openslide_read_associated_image(slidePtr,associatedImages{2});

% Display label image
figure
imshow(label(:,:,2:4))

% Close whole-slide image, note that the slidePtr must be removed manually
openslide_close(slidePtr)
clear slidePtr

% Unload library
openslide_unload_library