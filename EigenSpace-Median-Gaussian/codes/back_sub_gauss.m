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
%%
% Describe here your background subtraction method
u = 0;
sigma2 =1;
a = 0.03;
T =1.08;
se1 = strel('disk', 7);
se2 = strel('disk', 5);
%%%%%%%%%%%%%%%%%
TP = 0; FP = 0; FN = 0;
k = 0;


for i=795:NumImages
    
    I = ImSeq(:,:,i); 
    u = a.*I + (1-a)*u;
    sigma2 = (I-u).^2*a + (1- a).*sigma2;
    
    testIm = abs(I - u)./sqrt(sigma2);
    
    subplot(222)
    imshow(u, []);
    title('Background Model');
    Res = testIm > T;
    subplot(223)
    imshow(Res, [])
    title('Subtracted and thresholded');
    Res =medfilt2(Res, [5,5]);
    Res = bwareaopen(Res, 18);
    Res = imerode(Res, strel('disk', 3));
%     Res = imfill(Res, 'holes');
    Res = imdilate(Res, strel('square', 5));
    
 %% Bounding Box using regionprops
    s = regionprops(bwlabel(Res), 'centroid', 'BoundingBox');
    subplot(221)
    imshow(ImSeq(:,:,i), []);
    hold on;
    if ~isempty(s)
        for j=1:numel(s)
            rectangle('Position', s(j).BoundingBox, 'EdgeColor','g', 'LineWidth', 2);
            plot(s(j).Centroid(1), s(j).Centroid(2), 'r*');
        end
    end
    hold off
    title('tracking')

    
    subplot(224)
    imshow(Res);
    title('morphed image');
    pause(0.005)

    pause(0.003)
   
    TP = TP + sum(sum( Res & gndImSeq(:,:,i)));
    FP = FP + sum(sum((Res - gndImSeq(:,:,i)) > 0));
    FN = FN + sum(sum((gndImSeq(:,:,i) - Res) > 0));
end  

precision = TP/(TP + FP);
recall = TP/(TP + FN);

Fscore = 2*precision*recall/(precision + recall)




















