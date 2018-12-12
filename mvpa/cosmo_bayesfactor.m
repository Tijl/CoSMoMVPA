function ds_bf = cosmo_bayesfactor(ds,varargin)

    %

    %% deal with input
    defaults.h0mean = 0; %for decoding, 50% is chance
    defaults.prior = 'cauchy';
    defaults.cauchy_scale = 0.707;
    opt=cosmo_structjoin(defaults,varargin);
        
    %% subtract h0mean from ds
    samples = ds.samples-opt.h0mean;
    
    %% compute mean and se
    n = size(samples,1);
    MU = mean(samples,1);
    SE = std(samples,[],1)./sqrt(n);
    T = MU./SE;
    
    scale = opt.cauchy_scale;
    BF = bfun_t1smpbf_mat(T,n,scale);
    
    ds_bf = struct();
    ds_bf.samples = BF;
    ds_bf.sa.label = 'BF10';
    ds_bf.fa = ds.fa;
    ds_bf.a.fdim = ds.a.fdim;

function bf10 = bfun_t1smpbf_mat(t,n,scale)
%
% bf10 = t1smpbf(t,n,[r=0.707])
%
% Calculates JZS Bayes Factor for a one-sample t-test given t and sample size n.
% The optional input r is the scale factor which defaults to 0.707.
% This quantifies the evidence in favour of the alternative hypothesis. 
% See Rouder et al, 2009, Psychon Bull Rev for details.
% 
% original code obtained from https://doi.org/10.6084/m9.figshare.1357917.v1
% adapted to support a vector of t values
%

% Default scale factor
% if nargin < 3
%     scale = 0.707;
% end

% Function to be integrated
F = @(g,t,n,r) (1+n.*g.*r.^2).^(-1./2) .* (1 + t.^2./((1+n.*g.*r.^2).*(n-1))).^(-n./2) .* (2.*pi).^(-1./2) .* g.^(-3./2) .* exp(-1./(2.*g));

% Bayes factor calculation
bf01 = (1 + t.^2/(n-1)).^(-n/2) ./ integral(@(g) F(g,t,n,scale),0,Inf,'ArrayValued',true);

% Invert Bayes Factor
bf10 = 1 ./ bf01;
    