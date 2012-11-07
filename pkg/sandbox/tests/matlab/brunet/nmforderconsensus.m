function  [ordcons,clustid,ordindex,coph] = nmforderconsensus(consensus,kstart,kend)
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
%
% Reordering of consensus matrices using hierarchical clustering
%
% consensus : 3d array of consensus matrices 
%             dimensions : kend x M x M 
%
% kstart, kend : range of  consensus matrices to reorder
%
%
% ordcons  : 3d array of consensus matrices 
%            dimensions : kend x M x M
%            Values in ordcons(1:kstart-1,:,:) should be ignored
%
% clustid : dimensiosn kend x M
%           clustid(k,:) has cluster id for each sample (in original order)
%           Values in clustid(1:kstart-1,:) should be ignored
%
% ordindex : dimension kend x M 
%            ordindex (k,:) is the permutation vector 
%            for the kth consensus matrix 
%            Values in ordindex(1:kstart-1,:) should be ignored
%
% coph :  dimension kend
%         coph(k) has cophenetic coefficients for k 
%         Values in coph(1:kstart-1) should be ignored
%


[kmax,m,m]=size(consensus);

%check that kend is not out of range
if(kend>kmax) 
error('kend is out of range');
end

ordcons=zeros(kmax,m,m);
ordindex=zeros(m,kmax);
clustid=zeros(m,kmax);
coph=zeros(1,kmax);

for i=kstart:kend
u=reshape(consensus(i,:,:),m,m);
[ordcons(i,:,:),clustid(:,i),ordindex(:,i),coph(i)] = nmforderconsensus0(u,i);
end



function  [aordered,clust,ord,coph] = nmforderconsensus0(a,k)
%
% Jean-Philippe Brunet
% Cancer Genomics 
% 

[n,m]=size(a);
ordl=zeros(1,m);
incr=1;

uvec=a(1,2:end);

for i=2:n-1;
uvec=[uvec a(i,i+1:end)]; %get upper diagonal elements of consensus
end

y=1-uvec;                 % consensus are similarities, convert to distances
z=linkage(y,'average');   % use average linkage
coph=cophenet(z,y);

fig = figure('visible','off'); % turn off dendrogram plot
[h,t,ord]=dendrogram(z,0); % get permutation vector
close(fig)

clust=cluster(z,k);       % get cluster id 
aordered=a(ord,ord);
ord=ord';




