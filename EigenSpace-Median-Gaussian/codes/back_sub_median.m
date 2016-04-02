%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUAL TRACKING
% ----------------------
% Background Subtraction
% ----------------
% Date: september 2015
% Authors: Sandeep M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all


%%%%% LOAD THE IMAGES
%=======================

% Give image directory and extension
imPath = 'highway/input'; imExt = 'jpg';
groundPath = 'highway/groundtruth'; gExt = 'png';
% check if directory and files exist
if isdir(imPath) == 0
    error('USER ERROR : The image directory does not exist');
end

filearray = dir([imPath filesep '*.' imExt]); % get all files in the directory
groundarray = dir([groundPath filesep '*.' gExt]);
NumImages = size(filearray,1); % get the number of images
if NumImages < 0
    error('No image in the directory');
end

disp('Loading image files from the video sequence, please be patient...');
% Get image parameters
imgname = [imPath filesep filearray(1).name]; % get image name



I = imread(imgname); % read the 1st image and pick its size
VIDEO_WIDTH = size(I,2);
VIDEO_HEIGHT = size(I,1);

ImSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
gndImSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
for i=1:NumImages
    imgname = [imPath filesep filearray(i).name]; % get image name
    gndname = [groundPath filesep groundarray(i).name];
    ImSeq(:,:,i) = rgb2gray(imread(imgname)); % load image
    I = (imread(gndname));
    level = graythresh(I./255);
    gndImSeq(:,:,i) = im2bw(I, level); % load image
   
end
disp(' ... OK!');


% BACKGROUND SUBTRACTION
%=======================

% Describe here your background subtraction method

se = strel('square',12);
TP = 0; FP = 0;  FN = 0;
alpha = 0.05;
B = median(ImSeq(:,:, 1:470), 3);
for n =471:NumImages

    	
  
    diff = abs(ImSeq(:,:, n) - B);
    B = ImSeq(:,:, n).*alpha + B.*(1- alpha);
    Res = diff > 85;
    %% Morphological Operation
    Res = bwareaopen(Res, 10);
    Res = imfill(Res, 'holes');
    Res = imclose(Res, se);
    Res = imdilate(Res, se);
   %%%%%%%%%%%%%%%%%%%
    subplot(223)
    imshow(diff);
    title('difference image')
   %% Bounding Box using regionprops
    s = regionprops(bwlabel(Res), 'centroid', 'BoundingBox');
    subplot(221)
    imshow(ImSeq(:,:,n), []);
    hold on;
    if ~isempty(s)
        for i=1:numel(s)
            rectangle('Position', s(i).BoundingBox, 'EdgeColor','g', 'LineWidth', 2);
            plot(s(i).Centroid(1), s(i).Centroid(2), 'r*');
        end
    end
    hold off
    title('tracking')
    subplot(222)
    imshow( B, []);
    title('background image')
    
    subplot(224)
    imshow(Res);
    title('morphed image');
    pause(0.005)

   
    TP = TP + sum(sum( Res & gndImSeq(:,:,n)));
    FP = FP + sum(sum((Res - gndImSeq(:,:,n)) > 0));
    FN = FN + sum(sum((gndImSeq(:,:,n) - Res) > 0));
end  

precision = TP/(TP + FP);
recall = TP/(TP + FN);

Fscore = 2*precision*recall/(precision + recall)








