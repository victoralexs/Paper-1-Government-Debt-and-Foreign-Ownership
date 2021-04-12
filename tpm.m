cd '~/Google Drive/PhD Insper/Thesis/Paper 1 - Sovereign Default Holdings/Quant/Paper 1'

%%% TPM.m: Create a transition matrix for y shock %%%
%%%%%
%y_t = rho * y_t-1 + sigma_epsilon * epsilon_t, where epsilon_t follows a normal (0,1) distribution.  
clear all

rho = 0.70;
sigma_epsilon = 0.026;
epsilonvar = sigma_epsilon^2;

W = 4.0; %width of y_t grid  around its mean (0).
N1 = 200; % Grid points
T = 1e7; %Length of the time series of simulated

randn('state',0);%set seed
stdy= sigma_epsilon/sqrt(1-rho^2);
UB = W*stdy; 
ygrid = -UB: 2*UB / (N1-1) : UB;   % a grid for y spaced by 2*UB/N1-1 from -UB to UB, UB = W * stdy
PAI = zeros(N1); 

[xx,i0] = min(abs(ygrid));
y0 = ygrid(i0);

for t=2:T
y1 = y0*rho + randn*sigma_epsilon;
[xxx,i1] = min(abs(y1-ygrid));

PAI(i0,i1) = PAI(i0,i1) + 1;
y0 = y1;
i0 = i1;
end

%eliminate all rows and columns with all elements equal to zero
pai = PAI(sum(PAI,2)~=0,sum(PAI,1)~=0);
pai = pai ./ repmat(sum(pai,2),[1,size(pai,2)]);
keep = find(sum(PAI,2)~=0);
ygrid = ygrid(keep);
ygrid = ygrid(:);
N = length(ygrid); % N is the final lenght for ygrid, the grid for the exogenous process.

save tpm.mat  pai ygrid % Transition matrix