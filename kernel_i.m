function y=kernel_i(x,a,h,kernel)

%Integrals of Kernels: y = KK[dx K(x/h-a/h) ]
%
%====================   INPUT ============================================
%   x:          Vector (or matrix) of distances
%   a:          parameter
%   h:          Bandwidth
%   KERNEL:     Type of kernel funcion used in the estimate of the 
%                   sample constraints
%               Supported kernel functions 
%               Kernel choices: 
%               Epanechnikov:  'Epan'
%               Exponential:   'Expo'
%               Triweight:     'Triw'
%               Uniform:       'Unif'
%               Tricubic:      'Tric'
%               Bitriangular:  'Tria'
%               Spherical:    'Sphe'
%
%====================   OUTPUT ============================================
%
%   y:  Values of the integral of kernel function (vector or matrix, 
%   depending on x)
%
%==========================================================================

C=0; % Sta8era oloklhrwshs

if strcmp(kernel,'Expo')
    y = h * (x-a) .* (1 - exp(-abs(x-a)/h)) ./ abs(x-a) /2 + C;
    
elseif strcmp(kernel,'Gaus')
    y= h /2 * erf( (x-a)/h )  +C;
    
    
elseif strcmp(kernel,'Unif')
    y = x + C;
    y=y/2;
    
elseif strcmp(kernel,'Epan')
    y = x-a - (x-a).^3 / (3*h^2) + C;
    y= y*3/4;    
    
elseif strcmp(kernel,'Triw')
    y = x-a -(x-a).^7/(7*h^6)+(3*(x-a).^5)/(5*h^4)-(x-a).^3/h^2+C;
    y= 35*y/32;    
    
elseif strcmp(kernel,'Tric')
    y= -((x-a).^9.*abs(x-a))/(10*h^9)-(3*(x-a).^3.*abs(x-a))/(4*h^3)+(3.*(x-a).^7)/(7*h^6)+x-a +C;
    y= 1*y/1.1571;  
    
elseif strcmp(kernel,'Tria')
    y=-((x-a).*(3*h*abs(x-a)-x.^2 + 2*a*x - 3*h^2 - a^2))/(3*h^2);
    y=y*3/2;
    
elseif strcmp(kernel,'Sphe')
    y=((x-a).*((x-a).^2-6*h^2).*abs(x-a))/(8*h^3)+(x-a);
    y=y*4/3;
    
else
	sprintf('%s','Unrecognizable Kernel')
end







