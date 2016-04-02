%%
%Sandeep Manandhar and Anas MHANA
%Dec 30, 2015
%(dis)Similarity measure for direction models in GMM based DRA
%Background subtraction
%uB, France
%%
peakNum = 5;
forward_imPath = './../forwardRes'; imExt = 'jpg';
backward_imPath = './../backwardRes'; imExt = 'jpg';
%%equal number of samples in both models and same names of corresponding
%%frames assumed

if isdir(forward_imPath) == 0
    error('USER ERROR : The image directory does not exist');
end

filearray = dir([forward_imPath filesep '*.' imExt]); % get all files in the directory


NumImages = size(filearray,1); % get the number of images

if NumImages < 0
    error('No image in the directory');
end

disp('Loading image files from the video sequence, please be patient...');
% Get image parameters
imgname = [forward_imPath filesep filearray(3).name]; % get image name
I = imread(imgname);
ch = size(I,3) == 3;
VIDEO_WIDTH = size(I,2);
VIDEO_HEIGHT = size(I,1);
forwardSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
backwardSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
for i=1:NumImages
    f_imgname = [forward_imPath filesep filearray(i).name]; % get image name
    b_imgname = [backward_imPath filesep filearray(i).name];
    forwardSeq(:,:,i) = (imread(f_imgname));
    backwardSeq(:,:,i) = (imread(b_imgname));% load image
end

h = waitbar(0, '(dis)Similarity measure in progress...');
set(h, 'Name', 'GMM');
for i = 1:NumImages-3
    %     f_hist(:,:,i) =  0.5*(forwardSeq(:,:,i+2) - (forwardSeq(:,:,i+1) - forwardSeq(:,:,i))/2)+0.5*forwardSeq(:,:,i+1);
    %     b_hist(:,:,i) = 0.5*(backwardSeq(:,:,i+2) - (backwardSeq(:,:,i+1) - backwardSeq(:,:,i))/2)+0.5*backwardSeq(:,:,i+1);
    simi(i) = 4- (MI2(forwardSeq(:,:,i+2),forwardSeq(:,:,i+1)) + MI2(forwardSeq(:,:,i+2),forwardSeq(:,:,i)) + ...
        MI2(forwardSeq(:,:,i+3),forwardSeq(:,:,i+1)) + MI2(forwardSeq(:,:,i+3),forwardSeq(:,:,i))+...
              MI2(backwardSeq(:,:,i+2),backwardSeq(:,:,i+1)) + MI2(backwardSeq(:,:,i+2),backwardSeq(:,:,i)) + ...
              MI2(backwardSeq(:,:,i+3),backwardSeq(:,:,i+1)) + MI2(backwardSeq(:,:,i+3),backwardSeq(:,:,i)))/4;
             % MI2 function provided by http://www.flash.net/~strider2/matlab.htm


    waitbar(i/NumImages);
end
close(h);
%%
plot(simi, '-r', 'Linewidth', 2);
hold on;
plot(simi, 'ob', 'Markersize', 3);
plot(simi, 'ob', 'Markersize', 5);


[p, loc] = findpeaks(simi);
[~, i] = sort(p);
d = loc(i);
peakNum = 5;
plot(d(end-5:end), simi(d(end-5:end)), '*g', 'Markersize', 12);

