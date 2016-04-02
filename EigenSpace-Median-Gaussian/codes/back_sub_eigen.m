%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUAL TRACKING
% ----------------------
% Background Subtraction
% ----------------
% Date: september 2015
% Authors: SandeepM
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

%%
% BACKGROUND SUBTRACTION
%=======================

% % Describe here your background subtraction method
%%

se1 = strel('disk', 5);
se2 = strel('disk', 7);
%%%%%%%%%%%%%%%%%
%%
[r, c] = size(ImSeq(:,:,1));
stack = [];
for i=1:470
    stack = [stack  reshape(ImSeq(:,:,i), [],1)];
end
m = mean(stack, 2);
A = stack - repmat(m, 1, size(stack, 2));
%%
B = A';
COV =B*B';
[V, D] = eig(COV);
[D order] = sort(diag(D), 'descend');
d = inv(diag(D)).^(1/2);
V = V(:,order);
v = A*V*d;
%%
%[us ss vs] = svd(A, 0);
%%
%Uk = us(:, 1:2);
U = v(:, 1:2);
T = 50;
se = strel('line', 17,12);
seerode = strel('square', 12);
TP = 0; FP = 0; FN = 0;
for i = 470:NumImages
    Y = reshape(ImSeq(:,:,i), [], 1);
    p = U'*(Y - m);
    y = U*p + m;
    
    res =  abs(y - Y) > T;
    Res = reshape(res, r, c);
    subplot(223)
    imshow(Res, []);
   
    Res = bwareaopen(Res, 12);
    Res = imfill(Res, 'holes');
    Res = imclose(Res, se);
    s = regionprops(bwlabel(Res), 'centroid', 'BoundingBox');
    subplot(222);
    imshow(reshape(y, r, c), []);
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
    subplot(224)
    imshow(Res, []);
    pause(0.03)  
    
    TP = TP + sum(sum(Res & gndImSeq(:,:,i)));
    FP = FP + sum(sum((Res - gndImSeq(:,:,i)) > 0));
    FN = FN + sum(sum((gndImSeq(:,:,i) - Res) > 0));
end  

precision = TP/(TP + FP);
recall = TP/(TP + FN);

Fscore = 2*precision*recall/(precision + recall)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





















