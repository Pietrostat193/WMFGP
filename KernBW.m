function [ StpFun, CDF_Kern, zval ] = KernBW(data, kernel, Bandw, fig_Flag)
%	This function estimates the step function CDF for a timeseries and the 
%   non-parametric CDF approximation using Kernel functions.
%==========================================================================
% INPUTS:
%==========================================================================
%   data:       The sample values of the time series. N x 1
%   kernel:     The Kernel function to be used. 
%   Bandw:      The kernel bandwidth to be used
%   fig_Flag:   It is set to 0 if we don't want to generate a figure. 
%
%==========================================================================
%      *******  Kernel choices:         ******
%==========================================================================
%               Exponential:    'Expo'
%               Epanechnikov:   'Epan' 
%               Triweight:      'Triw'
%               Uniform:        'Unif'
%               Tricubic:       'Tric'
%               Bitriangular:   'Tria'
%               Sphehrical:     'Sphe'
%==========================================================================
%   OUTPUT
%==========================================================================
%   StpFun:     Step CDF function evaluated at 400  evenly placed steps.
%   CDF_Kern:   CDF estimated using the integration of the chosen kernel
%               function evaluated at 400  evenly placed steps.
%   zval:       Z values estimated at the 400 steps
%
%==========================================================================
%   FUNCTIONS used: kernel_s, kernel_i
%==========================================================================
% *************** EXAMPLE *************************************************
%==========================================================================
%   n=100; sz=[n,1]; 
%   pd = makedist('gam'); pd.a = 1; pd.b = 30;
%   rng(100); 
%   TS=random(pd,sz);
%   Kernelt='Sphe'; h=9;
%   [ StpFun, CDF_Kern, zval ] = KernBW_simple(TS, Kernelt, h, 1 );

Z=data; h=Bandw;

NT=400; %NT the 400 evaluation points
Zmin=(min(Z)); ZMAX=(max(Z));
x=linspace(Zmin,ZMAX,NT);

Za=sort(data); 
N=length(Za);

B(NT)=0;
for i=1:NT
    B(i)=sum(Z<=x(i)) / N;
end

if fig_Flag~=0
tt=figure;
axes1 = axes('Parent',tt);
plot(x,B,'-','Linewidth',1.5);
hold on
xlabel('Z')
ylabel('Cumulative probability')
set(axes1,'FontSize',14);
legend('step CDF','Location','southeast');
end


%% EXPO Or Gauss
% This is for infinitely extended kernels

if strcmpi(kernel,'Expo')==1 || strcmpi(kernel,'Gaus')==1

    h2=num2str(h);
    tit1=strcat(kernel,' h=',h2);

    t0=tic;
    B2 = zeros(NT, 1);
    for i=1:NT %Estimates the kernel in NT points from N points.
        Zsum=0;
        for j=1:N
            a=Za(j);
            ZZ=kernel_i(x(i),a,h,kernel) - kernel_i(x(1),a,h,kernel);
            Zsum=Zsum+ZZ;
        end
        B2(i)=sum( Zsum  / (N*h) );
    end
    t1=toc(t0);
    
    if fig_Flag~=0
        figure(tt)
        plot(x,B2/max(B2),'r');
        ylim([0 1.01]);
        legend('step CDF',tit1,'Location','southeast')
        set(gcf, 'Position',  [100, 100, 700, 500])
    end
end

%% compact

if strcmpi(kernel,'Expo')~=1 && strcmpi(kernel,'Gaus')~=1

    h2=num2str(h); tit1=strcat(kernel,' K, h=',h2);
    
    t0=tic;
    for i=1:NT % Estimate the kernel at NT points from N points.
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
    t1=toc(t0);

    if fig_Flag~=0    
        figure(tt)
        plot(x,B2,'r');
        hold on
        ylim([0 1.01]);

        legend('step CDF',tit1,'Location','southeast')
        set(gcf, 'Position',  [100, 100, 700, 500])
    end

end

CDF_Kern=B2/max(B2); zval=x; StpFun=B;


end

