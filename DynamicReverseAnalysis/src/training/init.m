%% Author : Sandeep Manandhar and Anas Mhanna
%%


%for each pixel sequence make gaussian mixture model
%get mean and covariance using EM
%Foreground pixel detection

%%


%%

clc
clear all
close all

%% read images
% imPath = 'car'; imExt = 'jpg';
%  imPath = 'LightSwitch'; imExt = 'bmp';
imPath = './backwardTrain'; imExt = 'bmp';

%%%%% LOAD THE IMAGES
%=======================
% check if directory and files exist
if isdir(imPath) == 0
    error('USER ERROR : The image directory does not exist');
end

filearray = dir([imPath filesep '*.' imExt]); % get all files in the directory
NumImages = size(filearray,1); % get the number of images
if NumImages < 0
    error('No image in the directory');
end

disp('Loading image files from the video sequence, please be patient...');
% Get image parameters
imgname = [imPath filesep filearray(1).name]; % get image name
I = imread(imgname);
ch = size(I,3) == 3;
VIDEO_WIDTH = size(I,2);
VIDEO_HEIGHT = size(I,1);

ImSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);

for i=1:NumImages
    imgname = [imPath filesep filearray(i).name]; % get image name
    
    if ch == 1
        ImSeq(:,:,i) = rgb2gray(imread(imgname)); % load image
    else
        ImSeq(:,:,i) = (imread(imgname));
    end
end
disp(' ... OK!');

%%
K = 2;
gaussianStat = zeros(VIDEO_HEIGHT, VIDEO_WIDTH,K*2);
tic;
m = reshape(ImSeq, [], size(ImSeq, 3));
%         [MU1 MU2 SIG1 SIG2] = EMGMM(m);
%

% [u1 u2 s1 s2 w1 w2] = EMGMM(m);
[u s w] = NGMM(m, 3);
%%

wt = w./s;
[sortedWt, idx] = sort(wt, 2);

for r = 1:size(wt, 1)
    sortedu(r,:) = u(r, idx(r,:));
    sorteds(r,:) = s(r, idx(r,:));
    sortedw(r,:) = w(r, idx(r,:));
end

U = reshape(sortedu, size(ImSeq,1), size(ImSeq, 2), 3);
S = reshape(sorteds, size(ImSeq,1), size(ImSeq, 2), 3);
W = reshape(sortedw, size(ImSeq,1), size(ImSeq, 2), 3);



toc

subplot(131); imshow(U(:,:,1), []);
title('low gaussian model')
subplot(132); imshow(U(:,:,2), []);
title('mid gaussian model')
subplot(133); imshow(U(:,:,3), []);
title('high gaussian model')
%%
%%file name for mean (u), variance(s), weights(W) followed by prexfix
% H for hybrid, F for forward and B for backward
save('uB.mat', 'U');
save('SB.mat', 'S');
save('WB.mat', 'W');
%%
