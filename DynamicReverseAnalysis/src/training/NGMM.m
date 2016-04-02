%% Author : Sandeep Manandhar and Anas Mhanna
%% Gaussian mixture model with 2 gaussians
%% 1-D
%% returns [mean1, mean2, sigma1, sigma2]
%% needs fix::try to remove NaN from iteration

function [u s w] = EMGMM(K,N)
% a = [ 134   136   134   134   132   133   132   132   135   131   131   131   134   132   130   134   133
%    132   132   132   128   125   130   128   130   131   131   130   128   129   131   131   130   133
%    129   132   132   133   132   131   132   132   130   132   135   126   131   132   129   132   132];
% K = a(:);
% K = K';
%%%
%% data synthesis
% w =zeros(size(K,1),N);
for i = 1:N
    w(:,i) = 1/N*ones(size(K,1),1);
    s(:,i) = std(K,0,2);
end

u = [mean(K,2) K(:,end) K(:,1)];
% 
% 
% w1 = 0.5*ones(size(K,1), 1);
% size(w1)
% w2 = 0.5*ones(size(K,1), 1);
% 
% u1 = mean(K, 2);
% size(u1)
% s1 = std(K,0,2);
% 
% s2 = s1;
% u2 = K(:,end);
% size(u2)





%%
h = waitbar(0, 'training...');
set(h, 'Name', 'Gaussian Mixture Mode');
for i = 1:600
    %% E-step
    %%weighted probabilities
    %% 

    % Y1 is the weighted liklihood 
    for j=1:N
      
         Y(:,:,j) = normpdfVec(K, u(:,j), s(:,j), w(:,j));
    end
    Ysum = sum(Y,3);
    for j = 1:N        
        t = Y(:,:,j);
        b(:,:,j) = t./Ysum;
    end
 

    %% M-step
    for j = 1:N
        t = b(:,:,j);
    w(:,j) = sum(t,2)/size(K,2);
   
    norml = (size(K,2)*w(:,j));
    
    up(:,j) = sum(t.*K,2)./norml;
    
    KminusMu1 = (K - repmat(up(:,j), 1, size(K,2)));
    
    
    sigN1 = sum(t.*KminusMu1.^2,2);
    
    
    sp(:,j) = sqrt(sigN1./norml);

    end
    if norm(up - u) <eps  
        disp('conv');
        break;
    end
    u = up;
    s = sp;
waitbar(i/200)
end
close(h);
% u1
%  s1
%  u2
%  s2
%  w1
%  w2
%  i
