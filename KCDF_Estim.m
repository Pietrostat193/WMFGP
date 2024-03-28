function [kcdf, nscores, h, Supplement] = KCDF_Estim(data, kernel)
%   KCDF_Estim generates a kernel-based estimate of the CDF of the
%   data using the BGK bandwidth. It also returns a vector of the normal
%   scores of the data values based on the kernel CDF. 
%==========================================================================
%           INPUT
%==========================================================================
%   data:           Vector of data values
%   kernel:         Name of kernel function (string variable)
%   Kernel choices: 
%       Epanechnikov:   'Epan'
%       Exponential:    'Expo'
%       Triweight:      'Triw'
%       Uniform:        'Unif'
%       Tricubic:       'Tric'
%       BiTriangular:   'Tria'
%       Spherical:      'Sphe'
%
%==========================================================================
%           OUTPUT
%==========================================================================
%   kcdf:           Vector of Kernel-Based CDF estimate (KCDF)
%   nscores:        Normal scores of the data values
%   h:              BGK bandwidth
%   Supplement:     Structure variable which contains the initial results
%                   of kernel-based CDF estimation at N > N(data) points
%   Supplement.pdf: KDE estimate of the PDF
%   Supplement.stairs: Staircase (empirical) estimate of the CDF
%   Supplement.zval: z-values where the KCDF is initially estimated
%   Supplement.nscores: normal scores of z values where the PDF
%               is estimated
%==========================================================================
%   FUNCTIONS used: kde, KernBW
%==========================================================================
%	REFERENCES
%   Bandwidth estimation using the Botev, Grotowski, Kroese rule.
%   Kernel density estimation via diffusion
%   Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
%   Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
%
%   Pavlides, A., Agou, V., & Hristopulos, D. T. (2021). Non-parametric 
%   Kernel-Based Estimation of Probability Distributions for Precipitation 
%   Modeling. arXiv preprint arXiv:2109.09961

ndisc = 2^14; % Number of uniformly distributed discretization points of 
%   the interval [zmin, zmax]

%   Minimum and maximum values for the interval on wihch the density
%   estimate is constructed.  Also, use 
%   >>help kde to see suggested values. 

zmin = eps; % For non-negative valued distributions
zmax = max(data) + 2; 

[h, pdf, ~, ~] = kde(data, ndisc, zmin, zmax); 
% h is the BGK bandwidth


%% CDF (Staircase and Kernel-derived) 

% Change 0 to 1 for figure generation
[ stairs_cdf, kcdf_estim, zval ] = KernBW(data, kernel, h, 0);

nscores1 = norminv(kcdf_estim); 

kcdf = interp1(zval, kcdf_estim, data); % CDF of the data points

nscores = norminv(kcdf); 
nscores(nscores==-Inf) = -3.5;
nscores(nscores==Inf) = 3.5;

Supplement.pdf = pdf; 
Supplement.stairs = stairs_cdf ; 
Supplement.zval = zval;
Supplement.nscores = nscores1;