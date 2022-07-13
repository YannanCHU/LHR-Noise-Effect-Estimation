clc;
clear;
close all;
addpath("..\..\myFunctions");

% ratio: 1:2200
% filename = "..\..\OS_Maps\Heathrow Airport\Colnbrook Poyle\map_7_1.png";
% filename = "..\..\OS_Maps\Heathrow Airport\Datchet\map_6_3.png";
% filename = "..\..\OS_Maps\Heathrow Airport\Harmondsworth_Harlington_Cranford\map_21_1.png";
% filename = "..\..\OS_Maps\Heathrow Airport\East Part1\map_8_16.png";
filename = "..\..\OS_Maps\Manchester Airport\Knutsford\map_1_4.png";

map = imread(filename); map = map(:,:,1:3);
worldFileName = getworldfilename(filename);
geoInfOfMap = worldfileread(worldFileName, 'planar', size(map));

% specify the threshold of the peak in 2D DFT results
peakFFTthreshold_low = 7000;
peakFFTthreshold_high = [];
% define the low frequency range
lowFreThreshold = 2048/50;
highFreThreshold = 2048/15;
% highFreThreshold = sqrt(sum(size(map).^2));
% specify the radius of circles in the mask
maskRadius = 7;
% the standard deviation of the 2D Gaussian kernel
gauSTD = 2;

mapGray = rgb2gray(map);
figure(1); imshow(mapGray); colormap gray; impixelinfo;

% 1 pixel = 0.25 meter  ==  2048 pixel = 512 meter
[mapWidth, mapLength, ~] = size(map);
resolution_x = geoInfOfMap.CellExtentInWorldX;
resolution_y = geoInfOfMap.CellExtentInWorldY;
alt = 11;

geoLatMin = geoInfOfMap.YWorldLimits(1);
geoLngMin = geoInfOfMap.XWorldLimits(1);
geoOrigin = [-geoLngMin, geoLatMin, alt];

% smapling spatial frequency
Fs_x = 1/resolution_x;      fx = linspace(-Fs_x/2,Fs_x/2,mapLength+1); fx = fx(1:end-1);
Fs_y = 1/resolution_y;      fy = linspace(-Fs_y/2,Fs_y/2,mapWidth+1); fy = fy(1:end-1);

mapBinary = (double(mapGray)>=223);
figure(2); imshow(mapBinary);  impixelinfo;

%% Apply the 2D DFT
mapBinaryFFT = fftshift(fft2(mapBinary));
mapBinaryMag = abs(mapBinaryFFT);
mapBinaryPhase = angle(mapBinaryFFT);

% use 2D Gaussian Filter to smoothen the DFT results
if isempty(gauSTD)
    mapBinaryMag2 = mapBinaryMag;
else
    mapBinaryMag2 = imgaussfilt(mapBinaryMag,gauSTD);
end
% figure(3); imagesc(log(mapBinaryMag2));  colormap('turbo'); impixelinfo; colorbar; caxis([6,9]); 

% remove low-frequency coefficients
cent = 1+([mapWidth, mapLength]/2);
m1 = 1:1:mapWidth;
n1 = 1:1:mapLength;
[x,y]=meshgrid(m1,n1);

freMat = sqrt((x-cent(2)).^2 + (y-cent(1)).^2);
indexOfFreWithinRange = (freMat>lowFreThreshold) & (freMat<highFreThreshold);

% plot the Magnitude spectrum and peaks
maxFFTValue = max(max(mapBinaryMag2));
[XMAX,IMAX,XMIN,IMIN] = extrema2(mapBinaryMag2.*indexOfFreWithinRange);
% eliminate peaks under minimum threshold
underThresh = (XMAX < peakFFTthreshold_low);

if isempty(peakFFTthreshold_high)
    beyondThresh = underThresh;
else
    aboveThresh = (XMAX > peakFFTthreshold_high);
    beyondThresh = underThresh | aboveThresh;
end
IMAX(beyondThresh) = [];
XMAX(beyondThresh) = [];

figure();
surf(mapBinaryMag2,'EdgeColor','none');         hold on;
zlim([0, 0.2*maxFFTValue]);
[y_index,x_index] = ind2sub(size(mapBinaryMag2),IMAX);
scatter3(x_index,y_index,XMAX,'r','filled');
title("The magnitude of 2D DFT results")
xlabel("Spatial Frequency in X");   ylabel("Spatial Frequency in Y");   zlabel("Magnitude");

figure();   subplot(1,2,1); 
imagesc(log(mapBinaryMag2)); colormap('turbo'); colorbar;    hold on;
scatter(x_index,y_index,'w','filled');
colormap("turbo");  colorbar; impixelinfo; hold off;    caxis([6,8]);
title("All peaks greater than FFT peak magnitude threshold");
daspect([1,1,1]);

% assume there are at most 3 valid candidate regions
if (length(XMAX) > 6)
    % three largest peak indexes
    [~,idx] = maxk(XMAX,6);
    XMAX = XMAX(idx);
    x_index = x_index(idx);
    y_index = y_index(idx);
end


subplot(1,2,2);  
imagesc(log(mapBinaryMag2)); title("log Magnitude of the 2D FFT of binary map");
colormap("turbo");  colorbar; impixelinfo; hold on;     caxis([6,8]);
plot(x_index(:), y_index(:), 'o', 'Color', 'w', 'MarkerFaceColor', 'w'); hold off;
daspect([1,1,1]);

figure();
subplot(1,2,1);
scatter3(x_index,mapWidth+2-y_index,XMAX,'r','filled'); hold on;
scatter3(cent(1), cent(2), 1, 'b', "filled"); hold off;
daspect([1,1,1]);                    title("Remaining 2D DFT peaks");
xlim([cent(1)-highFreThreshold,cent(1)+highFreThreshold]);     xlabel("X");
ylim([cent(1)-highFreThreshold,cent(1)+highFreThreshold]);     ylabel("Y");
zlabel("Magnitude");                 view(2);

subplot(1,2,2);
distanceVal = sqrt(sum(([x_index,mapWidth+2-y_index] - cent).^2,2));
scatter3(x_index, mapWidth+2-y_index, distanceVal, 'r', 'filled'); hold on;
scatter3(cent(1), cent(2), 1, 'b', "filled"); hold off;
daspect([1,1,1]);                    title("Distance to the centre point");
xlim([cent(2)-highFreThreshold,cent(2)+highFreThreshold]);     xlabel("X");
ylim([cent(2)-highFreThreshold,cent(2)+highFreThreshold]);     ylabel("Y");
zlabel("Distance");                  view(2);

%% construct the binary mask to retain peaks and their surrounding regions
region_valid_superimposed = zeros(mapWidth, mapLength);

numOfInvalidPeaks = 0;
for peakPairIndex = 1:ceil(length(XMAX)/2)
% for peakPairIndex = 3
    FFTMask = zeros(mapWidth, mapLength);
    % detect all harmonics
    peakIndex = 2*peakPairIndex;
    multipleVal_y = floor(abs(cent(1) / (y_index(peakIndex)-cent(1))));
    multipleVal_x = floor(abs(cent(2) / (x_index(peakIndex)-cent(2))));

    if multipleVal_y < multipleVal_x
        multipleVals = [(-multipleVal_y:1:-1), (1:1:multipleVal_y)];
    else
        multipleVals = [(-multipleVal_x:1:-1), (1:1:multipleVal_x)];
    end

    delta_x = x_index(peakIndex)-cent(2);
    delta_y = y_index(peakIndex)-cent(1);

    x_harm = cent(2) + multipleVals * delta_x;
    y_harm = cent(1) + multipleVals * delta_y;

    for j = 1:length(y_harm)
        distance_Sq = (y-y_harm(j)).^2 + (x-x_harm(j)).^2;
        FFTMask(distance_Sq <= maskRadius^2) = 1;
    end
    %     distance_Sq = (y-y_index(i)).^2 + (x-x_index(i)).^2;
    %     FFTMask(distance_Sq <= maskRadius^2) = 1;

    % figure(); imshow(FFTMask);
    
    mapBinaryMagFiltered = mapBinaryMag .* FFTMask;
    img_complex = mapBinaryMagFiltered .* exp(1j * mapBinaryPhase);
    img_ifft2 = abs(ifft2(ifftshift(img_complex)));
    figure(); subplot(1,3,1); imagesc(img_ifft2); title("Reconstructed Image"); colormap("turbo"); colorbar; daspect([1,1,1]); impixelinfo;
    
    img_ifft2_bw = img_ifft2 > max(img_ifft2(:)) * 0.4;
    subplot(1,3,2);  imshow(img_ifft2_bw); title("Reconstructed Binary Image"); daspect([1,1,1]);   impixelinfo;
    
    %% morphological operation to remove trivial regions
    img_ifft2_bw_morph = imdilate(img_ifft2_bw, strel('disk',11));
    % remove small white regions
    img_ifft2_filtered = smallRegionRemoval(img_ifft2_bw_morph, 600);
    % remove small black regions
    img_ifft2_filtered = smallRegionRemoval(imcomplement(img_ifft2_filtered), 3000);
    img_ifft2_filtered = imcomplement(img_ifft2_filtered);

    img_ifft2_bw_morph = imdilate(img_ifft2_filtered, strel('disk',5));
    region_candidate_BW = smallRegionRemoval(img_ifft2_bw_morph, 16000);
    region_candidate_BW = smallRegionRemoval(imcomplement(region_candidate_BW), 6000);
    region_candidate_BW = imcomplement(region_candidate_BW);
    region_candidate_lb = bwlabel(region_candidate_BW);
%     figure(); imshow(region_candidate_BW); title("After morphological operation"); daspect([1,1,1]);   impixelinfo;
    subplot(1,3,3);  imagesc(region_candidate_lb); daspect([1,1,1]); title("Labelled Binary Image"); impixelinfo;
    %% postprocessing - check the validity of detected regions
    info = regionprops(region_candidate_lb, 'Area', 'BoundingBox', 'Orientation', 'Centroid');
    % disp("peak value = " + XMAX(2*peakPairIndex) + "Area sum = " + sum([info.Area]));
    peakFFT_Mag_Area = 2*XMAX(2*peakPairIndex) / sum([info.Area]);

    spatialFrequency_1st = sqrt(sum( ([y_index(2*peakPairIndex), x_index(2*peakPairIndex)]-cent).^2 )) / size(mapBinary,1);
    spatialPeriod_1st = 1 / spatialFrequency_1st;
    orientation_1st = rad2deg(atan(-delta_y / delta_x));

    [region_candidate_lb_new, invalidFlag, MRs] = validityCheck(info, mapBinary, region_candidate_lb, spatialPeriod_1st, orientation_1st);
    
    if invalidFlag == false
        region_valid_BW = (region_candidate_lb_new>0);
        region_centroids = cat(1, info.Centroid);
        region_valid_centroids = region_centroids(MRs > 0,:);
        MRs_valid = {MRs(MRs>0)};
        % display the detected region
        figure(); tiledlayout(1,2,'TileSpacing','compact');
        titleStr = sprintf("Binary Mask\nMagnitude of fundamental component is %.2f,\n" + ...
            "Orientation of fundamental component is %.2f°", XMAX(2*peakPairIndex), rad2deg(atan(-delta_y / delta_x)));
        nexttile, imagesc(FFTMask); title(titleStr);  colormap('gray'); impixelinfo; daspect([1,1,1]);

        map_masked = labeloverlay(map, 255*imcomplement(region_valid_BW));
        titleStr = sprintf("Region(s) of interest\n(FFT magnitude / region area = %.4f)", peakFFT_Mag_Area);
        nexttile,    imshow(map_masked);     title(titleStr);
        text(region_valid_centroids(:,1)-20, region_valid_centroids(:,2), cell2mat(sprintfc("MR=%.4f", MRs_valid{:})))

        region_valid_superimposed = region_valid_superimposed+region_valid_BW;
    else
        numOfInvalidPeaks = numOfInvalidPeaks + 1;
    end
end
%%
if numOfInvalidPeaks == ceil(length(XMAX)/2)
    fprintf("No desired regions can be detected\n");
else
    map_masked2 = labeloverlay(map, 255*imcomplement((region_valid_superimposed>0)));
    figure();    imshow(map_masked2);     title("All regions of interest are superimposed"); impixelinfo;

    region_valid_superimposed_LB = bwlabel(region_valid_superimposed);
    centers = cat(1, regionprops(region_valid_superimposed_LB, 'Centroid').Centroid);

    geoCenters = {round([geoLngMin + 0.25*centers(:,1), geoLatMin + 0.25*(mapWidth-centers(:,2))])};
    text(centers(:,1), centers(:,2), cell2mat(sprintfc("%d,%d", geoCenters{:})));

    figure(13); imagesc(region_valid_superimposed_LB); impixelinfo;
    text(centers(:,1), centers(:,2), 'o');
end


% functions
function [region_candidate_lb, invalidFlag,MRs] = validityCheck(info, mapBinary, region_candidate_lb, spatialPeriod_1st, orientation_1st)
MRs = zeros(length(info),1);
invalidFlag = false;
numOfInvalidRegions = 0;

for peakIndex = 1:length(info)
% for peakIndex = 1
    boundInfo = info(peakIndex).BoundingBox;
    xrange = round(boundInfo(1)): 1: round(boundInfo(1)+boundInfo(3));
    yrange = round(boundInfo(2)): 1: round(boundInfo(2)+boundInfo(4));
    xrange(xrange < 1 | xrange > size(mapBinary,2)) = [];
    yrange(yrange < 1 | yrange > size(mapBinary,1)) = [];
    yrangeNum = numel(yrange);
    xrangeNum = numel(xrange);

    regionDetected = imcomplement(mapBinary(yrange, xrange));
    regionDetectedMag = abs(fftshift(fft2(regionDetected)));
    cent2 = size(regionDetectedMag)/2+1;
    
    m2 = 1:1:yrangeNum;
    n2 = 1:1:xrangeNum;
    [x2,y2]=meshgrid(n2,m2);
    freMat2 = sqrt(( (x2-cent2(2)) / xrangeNum ).^2 + ( (y2-cent2(1)) / yrangeNum ).^2);
    indexOfFreWithinRange2 = (freMat2 >= 1/50) & (freMat2 <= 1/10);

    [~, y_index2, x_index2] = maxkVals(regionDetectedMag.*indexOfFreWithinRange2, 2);
    spatialFrequency_2rd = sqrt(sum(((y_index2(2) - cent2(1)) / yrangeNum).^2 + ((x_index2(2) - cent2(2)) / xrangeNum).^2));
    spatialPeriod_2rd = 1 / spatialFrequency_2rd;

    delta_x2 = x_index2(2)-cent2(2);
    delta_y2 = y_index2(2)-cent2(1);
    orientation_2rd = rad2deg(atan(-delta_y2 / delta_x2));

    spatialPeriod_diff = abs(spatialPeriod_2rd - spatialPeriod_1st) / spatialPeriod_1st;
    orientation_diff = abs(mod(orientation_2rd,180) - mod(orientation_1st,180));
    [MR, ~] = maxRegularity(regionDetected);

%     figure();
%     imshow(regionDetected);
%     figure();
%     imagesc(log(regionDetectedMag)); colormap('turbo'); colorbar; impixelinfo; caxis([6,8]); daspect([1,1,1]); hold on;
%     scatter(x_index2, y_index2, 'w','filled'); hold off;
%     fprintf("The first measured spatial period is %2.3f, the second measured spatial period is %2.3f.\n", spatialPeriod_1st, spatialPeriod_2rd);
%     fprintf("The first measured orientation is %2.3f, the second measured orientation is %2.3f.\n", orientation_1st, orientation_2rd);
%     fprintf("The measured maximal regularity is %.4f", MR);
    
    if spatialPeriod_diff >= 0.05 || orientation_diff >= 10 || info(peakIndex).Area <= 50000 || MR < 0.015
        % set this invalid region as the background of the binary image
        region_candidate_lb(region_candidate_lb == peakIndex) = 0;
        numOfInvalidRegions = numOfInvalidRegions + 1;
        % fprintf("The difference is %2.3f, two measurements do not match. This detected region is invalid!\n", spatialPeriod_diff);
    else
        MRs(peakIndex) = MR;
%         fprintf("period_1st = %2.3f, period_2rd = %2.3f, orientation_1st = %2.3f, orientation_2rd = %2.3f\n", spatialPeriod_1st, spatialPeriod_2rd, orientation_1st, orientation_2rd);
%         fprintf("The spatial period difference is %2.3f, the orientation difference is %2.3f, two measurements do match. Valid.\n", spatialPeriod_diff, orientation_diff);
    end

end

if numOfInvalidRegions == length(info)
    invalidFlag = true;
end

end


function map_labelled = smallRegionRemoval(imageWithTrivialRegions, regionSizeThreshold)
% inputs:
% imageWithTrivialRegions: binary image with many small-size regions (white region, black ground)
% regionSizeThreshold: scalar which indicates the size threshold of small regions
% output:
% map_labelled: binary image where pixels in trivial regions are set as 0

% map_labelled = labelmatrix(bwconncomp(img_ifft2_bw_morph));
[map_labelled, ~] = bwlabel(imageWithTrivialRegions, 4);
map_region_Area = regionprops(map_labelled, 'Area');
small_region_idx = find([map_region_Area.Area] < regionSizeThreshold);

for i = 1:length(small_region_idx)
    [small_region_pixel_x, small_region_pixel_y] = find(map_labelled == small_region_idx(i));
    map_labelled(small_region_pixel_x, small_region_pixel_y) = 0;
end

map_labelled = logical(map_labelled);
end

function [coefMags, y_index, x_index] = maxkVals(inputMatrix, k)
% inputs:
% inputMatrix: the input matrix whose largest k values are found and
% located
% k: number of largest values
% outputs:
% coefMags: the found k largest values
% y_index: row indices of k largest values
% x_index: column indices of k largest values

[coefMags, indexs] = maxk(inputMatrix(:), k);
% index of rows
y_index = mod(indexs,size(inputMatrix,1));
% index of columns
x_index = ceil(indexs / size(inputMatrix,1));
end

