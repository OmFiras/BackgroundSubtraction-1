clear ImSeq;
imPath = './backwardTest'; imExt = 'bmp';





%% backward data
U = load('uB.mat');
U = U.U;

S = load('SB.mat');
S = S.S;

W = load('WB.mat');
W = W.W;

%%
%%%%% LOAD THE IMAGES
%=======================
% check if directory and files exist
if isdir(imPath) == 0
    error('USER ERROR : The image directory does not exist');
end

filearray = dir([imPath filesep '*.' imExt]); % get all files in the directory
% groundarray = dir([groundPath filesep '*.' gExt]);

NumImages = size(filearray,1); % get the number of images
% gndNumImages = size(groundarray,1);
if NumImages < 0
    error('No image in the directory');
end

disp('Loading image files from the video sequence, please be patient...');
% Get image parameters
imgname = [imPath filesep filearray(3).name]; % get image name
I = imread(imgname);
ch = size(I,3) == 3;
VIDEO_WIDTH = size(I,2);
VIDEO_HEIGHT = size(I,1);

testSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
% gndImSeq = zeros(VIDEO_HEIGHT, VIDEO_WIDTH, NumImages);
for i=1:NumImages
    imgname = [imPath filesep filearray(i).name]; % get image name
    
    if ch == 1
        testSeq(:,:,i) = rgb2gray(imread(imgname)); % load image
    else
        testSeq(:,:,i) = (imread(imgname));
    end
end


clc
close all
wt1 = W;  %ascending order of wt/st
M = U;
St = S;
se = strel('disk',2);
%%
xsd = 5;
alpha = 0.3;
FP =0; TP = 0; FN = 0;

T = 0.6;
GSz = 3;
se = strel('disk', 5);
for t = 1:NumImages
    
    
    tIm = testSeq(:,:,t);
    
    
    temp = zeros(size(wt1));
    threshmap = zeros(size(tIm));
    matchIdx = zeros(size(wt1));
    
    %start searching from high rankin gaussian
    %1 - only high ranking gaussian involved
    %2 - medium and high ranking gaussians involved
    
    
    for k = GSz:-1:1 %%dont consider the last one[sum always equal to 1 if added]
        temp(:,:,k) = sum(wt1(:,:,GSz:-1:k),3);
        threshmap = max(threshmap, (temp(:,:,k) > T)*(k) );
    end
    
    %%how many gaussians map = threshmap
    threshmap = (GSz +1 - threshmap);
    threshmap(threshmap>GSz) = 0;
    
    for k = GSz:-1:1
        matchIdx(:,:,k) = (abs(tIm - M(:,:,k)) <= xsd*St(:,:,k));
    end
    
    
    firstmatch = matchIdx;
    firstmatch(:,:,2) = matchIdx(:,:,2) & ~matchIdx(:,:,3); %didnt match with 3 but with matched with 2
    firstmatch(:,:,1) = matchIdx(:,:,1) & (~matchIdx(:,:,2) & ~matchIdx(:,:,3));
    unmatch = ~matchIdx(:,:,1) & (~matchIdx(:,:,2) & ~matchIdx(:,:,3));
    
    wt1 = wt1 + alpha*(firstmatch - wt1);
    
    
    
    M(:,:,1) = unmatch.*(tIm) +  ~unmatch.*M(:,:,1);    %M(:,:,1) = meanImg;
    St(:,:,1) = unmatch.*(10) +  ~unmatch.*St(:,:,1);    %M(:,:,1) = meanImg;
    wt1(:,:,1) = unmatch.*(0.2) + ~unmatch.*wt1(:,:,1);
    

    sumWt = sum(wt1,3);
    
    wt1= bsxfun(@rdivide, wt1, sumWt);
    
    %the gaussian that got matched will be updated
    
    
    for nn = 1:2
        rho(:,:,nn) = alpha*exp(-0.5*((tIm - M(:,:,nn))./St(:,:,nn)).^2)./(sqrt(2*pi)*St(:,:,nn));
        rho(:,:,nn) = rho(:,:,nn).*firstmatch(:,:,nn);
      
        rho((rho >1)) = 1;
        rho(find(isnan(rho))) = 0;
        
        M(:,:,nn) = (1 - rho(:,:,nn)).*M(:,:,nn) + rho(:,:,nn).*tIm;
        St(:,:,nn) = sqrt((1 - rho(:,:,nn)).*St(:,:,nn).^2 +...
            rho(:,:,nn).*(tIm - M(:,:,nn)).^2);
    end
    
    rank = wt1./St;
    
    rank = reshape(rank, [], 3);
    wt1  = reshape(wt1, [], 3);
    M    = reshape(M, [], 3);
    St   = reshape(St, [], 3);
    
    [sd, idx] = sort(rank, 2);
    
    for r = 1:size(rank, 1)
        sortedu(r,:) = M(r, idx(r,:));
        sorteds(r,:) = St(r, idx(r,:));
        sortedw(r,:) = wt1(r, idx(r,:));
    end
    
    M = reshape(sortedu, size(testSeq,1), size(testSeq, 2), 3);
    St = reshape(sorteds, size(testSeq,1), size(testSeq, 2), 3);
    wt1 = reshape(sortedw, size(testSeq,1), size(testSeq, 2), 3);

    matching = zeros(size(matchIdx));
    for k = GSz:-1:1
        matching(:,:,k) =( (threshmap==(GSz-k+1)) & (firstmatch(:,:,k) ));
    end
    %new pixles that did not match with any models
    fg = ~matching(:,:,1) & (~matching(:,:,2) & ~matching(:,:,3)); 

    
    subplot(231);
    imshow(M(:,:,1),[]);
    subplot(232);
    imshow(M(:,:,2),[]);
    subplot(233);
    imshow(M(:,:,3),[]);
    
        
    subplot(234)
    imshow(tIm,[]);
        
    subplot(235)
    fgo = bwareaopen(medfilt2(fg), 12);

    storeIm = mat2gray(fgo.*tIm);
    imwrite(storeIm, strcat('./../backwardRes/b',num2str(t), '.jpg'));
    imshow(storeIm,[]);
    title('foreground')

    subplot(236)
    imshow(matching(:,:,3),[]); %%background is wh

    
    
    
    pause(0.04)

end
















