% M. H. Van Benthem and M. R. Keenan, J. Chemometrics 2004; 18: 441-450
%
% Given A and C this algorithm solves for the optimal 
% K in a least squares sense, using that
%      A = C*K 
% in the problem
%      min ||A-C*K||, s.t. K>=0, for given A and C.
%
function [K, Pset] = fcnnls(C, A)
% NNLS using normal equations and the fast combinatorial strategy
%
% I/O: [K, Pset] = fcnnls(C, A);
% K = fcnnls(C, A);
%
% C is the nObs x lVar coefficient matrix
% A is the nObs x pRHS matrix of observations
% K is the lVar x pRHS solution matrix
% Pset is the lVar x pRHS passive set logical array
%
% M. H. Van Benthem and M. R. Keenan
% Sandia National Laboratories
%
% Pset: set of passive sets, one for each column
% Fset: set of column indices for solutions that have not yet converged
% Hset: set of column indices for currently infeasible solutions
% Jset: working set of column indices for currently optimal solutions
%
% Check the input arguments for consistency and initialize
error(nargchk(2,2,nargin))
[nObs, lVar] = size(C);
if size(A,1)~= nObs, error('C and A have imcompatible sizes'), end
pRHS = size(A,2);
W = zeros(lVar, pRHS);
iter=0; maxiter=3*lVar;
% Precompute parts of pseudoinverse
CtC = C'*C; CtA = C'*A;
% Obtain the initial feasible solution and corresponding passive set
K = cssls(CtC, CtA);
Pset = K > 0;
K(~Pset) = 0;
D = K;
Fset = find(~all(Pset));
% Active set algorithm for NNLS main loop
oitr=0; % HKim
while ~isempty(Fset)
    
    oitr=oitr+1; if oitr > 5, fprintf('%d ',oitr);, end % HKim
    
    % Solve for the passive variables (uses subroutine below)
    K(:,Fset) = cssls(CtC, CtA(:,Fset), Pset(:,Fset));
    % Find any infeasible solutions
    Hset = Fset(find(any(K(:,Fset) < 0)));
    % Make infeasible solutions feasible (standard NNLS inner loop)
    if ~isempty(Hset)
      nHset = length(Hset);
      alpha = zeros(lVar, nHset);
      while ~isempty(Hset) & (iter < maxiter)
            iter = iter + 1; 
            alpha(:,1:nHset) = Inf;
            % Find indices of negative variables in passive set
            [i, j] = find(Pset(:,Hset) & (K(:,Hset) < 0));
            if isempty(i), break, end
            hIdx = sub2ind([lVar nHset], i, j);
            if nHset==1, % HKim
                negIdx = sub2ind(size(K), i, Hset*ones(length(j),1)); %HKim
            else % HKim
               negIdx = sub2ind(size(K), i, Hset(j)');
            end % HKim
            alpha(hIdx) = D(negIdx)./(D(negIdx) - K(negIdx));
            [alphaMin,minIdx] = min(alpha(:,1:nHset));
            alpha(:,1:nHset) = repmat(alphaMin, lVar, 1);
            D(:,Hset) = D(:,Hset)-alpha(:,1:nHset).*(D(:,Hset)-K(:,Hset));
            idx2zero = sub2ind(size(D), minIdx, Hset);
            D(idx2zero) = 0;
            Pset(idx2zero) = 0;
            K(:, Hset) = cssls(CtC, CtA(:,Hset), Pset(:,Hset));
            Hset = find(any(K < 0)); nHset = length(Hset);
      end
   end%if
   % Make sure the solution has converged
   %if iter == maxiter, error('Maximum number iterations exceeded'), end
   % Check solutions for optimality
   W(:,Fset) = CtA(:,Fset)-CtC*K(:,Fset);
   Jset = find(all(~Pset(:,Fset).*W(:,Fset) <= 0));
   Fset = setdiff(Fset, Fset(Jset));
   % For non-optimal solutions, add the appropriate variable to Pset
   if ~isempty(Fset)
       [mx, mxidx] = max(~Pset(:,Fset).*W(:,Fset));
       Pset(sub2ind([lVar pRHS], mxidx, Fset)) = 1;
       D(:,Fset) = K(:,Fset);
   end
end
end

% ****************************** Subroutine****************************
function [K] = cssls(CtC, CtA, Pset)
% Solve the set of equations CtA = CtC*K for the variables in set Pset
% using the fast combinatorial approach
K = zeros(size(CtA));
if (nargin == 2) || isempty(Pset) || all(Pset(:))
    K = CtC\CtA;
    %K=pinv(CtC)*CtA;
else
   [lVar pRHS] = size(Pset);
   codedPset = 2.^(lVar-1:-1:0)*Pset;
   [sortedPset, sortedEset] = sort(codedPset);
   breaks = diff(sortedPset);
   if size(breaks)(1) != 0
   breakIdx = [0 find(breaks) pRHS];
   else
   breakIdx = [0 pRHS];
   end
   for k = 1:length(breakIdx)-1
     cols2solve = sortedEset(breakIdx(k)+1:breakIdx(k+1));
     vars = Pset(:,sortedEset(breakIdx(k)+1));
     K(vars,cols2solve) = CtC(vars,vars)\CtA(vars,cols2solve);
     %K(vars,cols2solve) = pinv(CtC(vars,vars))*CtA(vars,cols2solve);
  end
end
