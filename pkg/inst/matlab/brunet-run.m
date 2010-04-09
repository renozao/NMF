%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab/Octave code to run NMF using Brunet et al.
% original algorithm (in file brunet.m).
%
% The original function was adapted to accept arguments to set the initial
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
[W,H] = nmf(target, rank, true, W0, H0);

% save result
save 'ref.brunet.oct' W H;
