function [ lookup] = Gen_Lookup( data, h, kernel, Npts )
%   [ lookup] = Gen_Lookup( data, h, kernel )
%   This function generates a lookup Table with 4000 points used to invert 
%   the normal scores transform 
%
%==========================================================================
%           INPUT
%==========================================================================

%	data:   The timeseries on which the simulation is based on
%   h:      Kernel Bandwidth
%	kernel: Kernel Model 
%   Supported kernel functions:
%     Epanechnikov:   'Epan'
%     Triweight:      'Triw'
%     Uniform:        'Unif'
%     Tricubic:       'Tric'
%     BiTriangular:   'Tria'
%     Sphehrical:     'Sphe'
%	Npts: (optional) The number of  points used in lookup table
%
%==========================================================================
%           OUTPUT
%==========================================================================
%	lookup: Look Up Table for the CDF derived from data
%==========================================================================
%   FUNCTIONS used: kernel_i
%==========================================================================
if nargin ==3
    NN=4000;
elseif nargin ==4
    NN = Npts;
else
    error('Wrong number of input arguments (3-4 needed)');
end
N=length(data); 
Za=sort(data);
x=linspace(min(data)-h, max(data)+h, NN); % Vector of values where the CDF is estimated

%% Step 1: Estimate CDF at all NN points

B2 = zeros(NN, 1);
for i=1:NN %Estimates the kernel in NN points from N sample points.
    Zsum=0;
    for j=1:N
        a=Za(j);
        if x(i)<a-h
            ZZ=0;
        elseif x(i)>a+h
            ZZ=kernel_i(a+h,a,h,kernel) - kernel_i(a-h,a,h,kernel);
        else
            ZZ=kernel_i(x(i),a,h,kernel) - kernel_i(a-h,a,h,kernel);
        end
        Zsum=Zsum+ZZ;
    end
    B2(i)=sum( Zsum  / (N*h) );
end

%% Step 2: Look Up Table

% k=sum(x'<0);
% B2(1:k)=[];
% x(1:k)=[]; 
% x(1)=eps;
kcdfx = B2;
lookup=[x(:) kcdfx(:)];

end

