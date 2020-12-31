function angle=Laplace_distribution(dev_standard)
% This function generates a Laplacian random variable with zero mean and
% standard deviation dev_standard.


%% INPUT PARAMETER

% dev_standard: standard deviation of Laplacian random variable in radiants

%% OUTPUT PARAMETER

% angle: Laplacian random variable with zero mean and standard deviation
% dev_standard

c=1/((sqrt(2)*dev_standard)*(1-exp(-sqrt(2)*pi/dev_standard)));
d=sqrt(2)/dev_standard;

x=rand; % uniform random variable in [0,1]

if x<=1/2
    angle=1/d*log(d/c*x+exp(-d*pi));
else
    angle=-1/d*log(1-(x-1/2)*d/c);
end
end
