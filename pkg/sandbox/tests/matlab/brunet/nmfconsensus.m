function consensus = nmfconsensus(a,kstart,kend,nloop,verbose)
%
% Jean-Philippe Brunet
% Cancer Genomics 
% The Broad Institute
% brunet@broad.mit.edu
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
% or functionality. 
%
% Model selection for NMF
%
% a (n,m) : N (genes) x M (samples) original matrix 
%           numerical data only. Must be positive.
% 
% kstart, kend : range of values of k to test consensus for.
%
% nloop : number of initial conditions per k 
%         (start with 10 or 20, check with more)
%
% verbose : prints iteration count and changes in connectivity matrix elements
%           if not set to 0 
%
% consensus : 3d array of consensus matrices 
%             dimensions : kend x M x M 
%             Values in consensus(1:kstart-1,:,:) should be ignored
%              

% test for negative values in v
if min(min(a)) < 0
error('matrix entries can not be negative');
return
end
if min(sum(a,2)) == 0
error('not all entries in a row can be zero');
return
end

[n,m]=size(a);

consensus=zeros(kend,m,m);
conn=zeros(m,m); 
incr=0;

for j=kstart:kend

if verbose fprintf(1,'rank %d\n',j), end

connac=zeros(m,m); 

for iloop=1:nloop;

[target, w, h] = load_test_data('data.5000.3.38', iloop)
%file = sprintf("data/W.%i.txt",iloop);
%w = load(file);
%file = sprintf("data/H.%i.txt",iloop);
%h = load(file);

incr=incr+1;

if verbose fprintf(1,' iteration %d\n',iloop), end 

[w,h]=nmf(a,w,h,j,verbose);

%
% compute consensus matrix
%
conn=nmfconnectivity(h); 
connac=connac+conn; % accumulate connectivity matrices

end

consensus(j, :, :)=connac/nloop; %average

end
end



function conn= nmfconnectivity(h)
%
% Jean-Philippe Brunet
% Cancer Genomics 6/10/03
%
mm=size(h);
k=mm(1);
m=mm(2);


% compute m x m matrix which is 1 if samples are together, 0 elsewhere

% determine sample assignment by its largest metagene expresion value
[y,index]=max(h,[],1); 

mat1=repmat(index,m,1); % spread index down
mat2=repmat(index',1,m); % spread index right

conn=mat1==mat2; % 1 when for pair of samples with same assignement



