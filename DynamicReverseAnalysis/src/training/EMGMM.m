%% Author : Sandeep Manandhar and Anas Mhanna
%% Gaussian mixture model with 2 gaussians
%% 1-D
%% returns [mean1, mean2, sigma1, sigma2]
%% needs fix::try to remove NaN from iteration

function [u1 u2 s1 s2 w1 w2] = EMGMM(K)
% a = [ 134   136   134   134   132   133   132   132   135   131   131   131   134   132   130   134   133
%    132   132   132   128   125   130   128   130   131   131   130   128   129   131   131   130   133
%    129   132   132   133   132   131   132   132   130   132   135   126   131   132   129   132   132];
% K = a(:);
% K = K';
%%%
%% data synthesis


w1 = 0.5*ones(size(K,1), 1);
size(w1)
w2 = 0.5*ones(size(K,1), 1);

u1 = mean(K, 2);
size(u1)
s1 = std(K,0,2);

s2 = s1;
u2 = K(:,end);
size(u2)






%%
for i = 1:400
    %% E-step
    %%weighted probabilities
  tic
    
    %% 
    
    % Y1 is the weighted liklihood 
    
    Y1 = normpdfVec(K, u1, s1, w1);
    Y2 = normpdfVec(K, u2, s2, w2);
    
    b = Y1./(Y1 + Y2) ;
    a = 1 - b;

    %% M-step
    w1 = sum(b,2)/size(K,2); %%sum along the col
    w2 = sum(a,2)/size(K,2);
    norml = (size(K,2)*w1);
    u1p = sum(b.*K,2)./norml;
    u2p = sum(a.*K,2)./norml;
    
    KminusMu1 = (K - repmat(u1p, 1, size(K,2)));
    KminusMu2 = (K -repmat(u2p, 1, size(K,2)));
    
    sigN1 = sum(b.*KminusMu1.^2,2);
    sigN2 = sum(a.*KminusMu2.^2,2);
    
    s1p = sqrt(sigN1./norml);
    s2p = sqrt(sigN2./norml);
    
    if norm(u1p - u1) <eps  
        disp('conv');
        break;
    end
    
    
  
    u1 = u1p;
    u2 = u2p;
    s1 = s1p;
    s2 = s2p;
toc
end

% u1
%  s1
%  u2
%  s2
%  w1
%  w2
%  i
