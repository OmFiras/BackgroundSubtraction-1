function y = normpdfVec(x,mu,sigma,w)



if nargin<1
    error(message('stats:normpdf:TooFewInputs'));
end
if nargin < 2
    mu = zeros(size(x,1), 1);
end
if nargin < 3
    sigma = ones(size(x,1), 1);
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;
MU = repmat(mu, 1, size(x, 2));
SIG = repmat(sigma, 1, size(x,2));
Wt = repmat(w, 1, size(x,2));

% try
    y = Wt.*exp(-0.5 * ((x-MU)./SIG).^2)./(sqrt(2*pi)*SIG);
    y(find(isnan(y)==1)) = 1;
% catch
%     error(message('stats:normpdf:InputSizeMismatch'));
end

