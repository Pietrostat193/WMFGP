function [zestim] = Kernel_invNS(normscores, lookup )
%==========================================================================
%           INPUT
%==========================================================================
%   normscores:     Vector of normal scores
%   lookup:         N x 2 array containing z values and respective CDF 
%                   estimates at a number of points (default 10000) 
%==========================================================================
%           OUTPUT
%==========================================================================
%   zestim:         Vector of z values corresponding to normscores
%==========================================================================
%==========================================================================
%   FUNCTIONS used: interp1
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

cdf_pred = normcdf(normscores); 

zvalues = lookup(:, 1); 

cdfz = lookup(:, 2); 

IDX = knnsearch(cdfz(:), cdf_pred(:));  %IDX the corresponding index 
zestim = zvalues(IDX); 

%   Sort the CDF values corresponding to predictions
%[cdfpred_sort, isort] = sort(cdf_pred);
% zestim = interp1(cdfz, zvalues, cdfpred_sort);
%   Undo the sorting
%w(isort) = zestim; 
%zestim = w; 

