%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB/Octave code to run NMF using Brunet et al. original algorithm 
% (see in file brunet.m).
% The objective is to be able to compare the results/performances from the 
% NMF package and the original MATLAB code.
%
% The original function 'nmf' was adapted to accept arguments to set the initial
% values for the matrix factors W and H.
%
% Original MATLAB codes can be found at:
% http://www.broadinstitute.org/mpr/publications/projects/NMF/nmf.m
% http://www.broadinstitute.org/publications/broad872
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean workspace
clear -all;

% load Brunet et al. algorithm -> defines the 'nmf' function
source('brunet.m');

% load test data
target = load('target.txt');
W0 = load('W.txt');
H0 = load('H.txt');
[n, rank] = size(W0);

% run algorithm
[total, user, sys] = cputime();
[W,H] = nmf(target, rank, false, W0, H0);
[total2, user2, sys2] = cputime();
elapsed = total2 - total;
user = user2 - user;
sys = sys2 - sys;
 
% save result
save 'ref.brunet.oct' W H elapsed user sys;
