function res =CSN_BM(Mu,Sigma,Gamma,Nu,Delta,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [res, t_wi, t_inv] =CSN_BM(Mu,Sigma,Gamma,Nu,Delta,num)
%   
%   Mu,Sigma,Gamma,Nu,Delta: parameters in CSN 
%   num: number of realisations
%
%
%   res: returned realisations
%
%
%   notes:   see dahoiv.net/master
%
% Created by Daniel Hoyer Iversen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % tic
 t=cputime;
  n2=length(Mu);
  n1=length(Nu);
  Gamma_vt=Gamma*Sigma;
  Sigma_v=Delta+Gamma*Sigma*Gamma';
  GSinv=Gamma_vt'*inv(Sigma_v);
  
  A1=chol(Sigma_v,'lower');
  A2=chol(Sigma-GSinv*Gamma_vt,'lower');
  
  res=zeros(num,n2);
  i=0;
  k=0;
  while i < num,
    i=i+1;
    z1=randn(n1,1);
    v=A1*z1-Nu;
    while any(v<0)
      k=k+1;
      z1=randn(n1,1);
      v=A1*z1-Nu;
    end
	cputime-t;
    z2=randn(n2,1);
    res(i,:)=(Mu+GSinv*(v+Nu)+A2*z2)';
  end
%  k/num
%  mean(res)
%  cputime-t
%  toc
end
  
