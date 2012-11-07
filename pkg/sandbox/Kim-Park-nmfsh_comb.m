%
% SNMF/R  
%
% Author: Hyunsoo Kim and Haesun Park, Georgia Insitute of Technology
%
% Reference: 
%
%   Sparse Non-negative Matrix Factorizations via Alternating 
%   Non-negativity-constrained Least Squares for Microarray Data Analysis
%   Hyunsoo Kim and Haesun Park, Bioinformatics, 2007, to appear.
%
% This software requires fcnnls.m, which can be obtained from 
% M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
%
% NMF: min_{W,H} (1/2) || A - WH ||_F^2 s.t. W>=0, H>=0 
% SNMF/R: NMF with additional sparsity constraints on H
%
%   min_{W,H} (1/2) (|| A - WH ||_F^2 + eta ||W||_F^2 
%                + beta (sum_(j=1)^n ||H(:,j)||_1^2))
%                s.t. W>=0, H>=0 
%
% A: m x n data matrix (m: features, n: data points)
% W: m x k basis matrix
% H: k x n coefficient matrix
%
% function [W,H,i]=nmfsh_comb(A,k,param,verbose,bi_conv,eps_conv)
%
% input parameters:
%   A: m x n data matrix (m: features, n: data points)
%   k: desired positive integer k
%   param=[eta beta]:  
%      eta (for supressing ||W||_F)
%         if eta < 0, software uses maxmum value in A as eta. 
%      beta (for sparsity control)
%         Larger beta generates higher sparseness on H.
%         Too large beta is not recommended. 
%   verbos: verbose = 0 for silence mode, otherwise print output
%   eps_conv: KKT convergence test (default eps_conv = 1e-4)
%   bi_conv=[wminchange iconv] biclustering convergence test 
%        wminchange: the minimal allowance of the change of 
%        row-clusters  (default wminchange=0)
%        iconv: decide convergence if row-clusters (within wminchange)
%        and column-clusters have not changed for iconv convergence 
%        checks. (default iconv=10)
%
% output:
%   W: m x k basis matrix
%   H: k x n coefficient matrix
%   i: the number of iterations
%
% sample usage:
%  [W,H]=nmfsh_comb(amlall,3,[-1 0.01],1);
%  [W,H]=nmfsh_comb(amlall,3,[-1 0.01],1,[3 10]); 
%     -- in the convergence check, the change of row-clusters to
%        at most three rows is allowed.
%
%
function [W,H,i]=nmfsh_comb(A,k,param,verbose,bi_conv,eps_conv)

if nargin<6, eps_conv=1e-4;, end
if nargin<5, bi_conv=[0 10];, end
if nargin<4, verbose=0;, end
if nargin<3, error('too small number of input arguments.');, end
maxiter = 20000; % maximum number of iterations

[m,n]=size(A); erravg1=[];   
eta=param(1); beta=param(2); 
maxA=max(A(:)); if eta<0, eta=maxA; end
eta2=eta^2;
wminchange=bi_conv(1); iconv=bi_conv(2);
if verbose, fprintf('SNMF/R k=%d eta=%.4e beta (for sparse H)=%.4e wminchange=%d iconv=%d\n',...
        k,eta,beta,wminchange,iconv);, end

idxWold=zeros(m,1); idxHold=zeros(1,n); inc=0; 

% initialize random W
W=rand(m,k);
W=W./repmat(sqrt(sum(W.^2,1)),m,1);  % normalize 

I_k=eta*eye(k); betavec=sqrt(beta)*ones(1,k); nrestart=0;
for i=1:maxiter

  % min_h ||[[W; 1 ... 1]*H  - [A; 0 ... 0]||, s.t. H>=0, for given A and W.
  H = fcnnls([W; betavec],[A; zeros(1,n)]);

  if find(sum(H,2)==0), 
      fprintf('iter%d: 0 row in H eta=%.4e restart!\n',i,eta);
      nrestart=nrestart+1;
      if nrestart >= 10, 
          fprintf('[*Warning*] too many restarts due to too big beta value...\n');, 
          break;
      end
      idxWold=zeros(m,1); idxHold=zeros(1,n); inc=0; 
      W=rand(m,k); W=W./repmat(sqrt(sum(W.^2,1)),m,1);  % normalize 
      continue;
  end
  
  % min_w ||[H'; I_k]*W' - [A'; 0]||, s.t. W>=0, for given A and H.
  Wt=fcnnls([H'; I_k],[A'; zeros(k,m)]); W=Wt';   

  % test convergence every 5 iterations
  if(mod(i,5)==0)  | (i==1)
      [y,idxW]=max(W,[],2);  [y,idxH]=max(H,[],1);    
      changedW=length(find(idxW ~= idxWold)); changedH=length(find(idxH ~= idxHold));
      if (changedW<=wminchange) & (changedH==0), inc=inc+1;, else inc=0;, end
      
      resmat=min(H,(W'*W)*H-W'*A+beta*ones(k,k)*H); resvec=resmat(:);
      resmat=min(W,W*(H*H')-A*H'+eta2*W); resvec=[resvec; resmat(:)];      
      conv=norm(resvec,1); %L1-norm      
      convnum=length(find(abs(resvec)>0));
      erravg=conv/convnum;
      if i==1, erravg1=erravg;, end
      
      if verbose | (mod(i,1000)==0)  % prints number of changing elements 
         fprintf('\t%d\t%d\t%d %d --- erravg1: %.4e erravg: %.4e\n',...
             i,inc,changedW,changedH,erravg1,erravg);
      end
      if (inc>=iconv) & (erravg<=eps_conv*erravg1), break, end 
      idxWold=idxW; idxHold=idxH; 
  end
  
end

return;