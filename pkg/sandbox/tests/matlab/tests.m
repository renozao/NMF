%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab/Octave code to perform tests for package NMF in R
%
% - tests to compare with already published algorithms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read matrices from files
function [target, W0, H0]=load_test_data(dir, n)
	target = load(sprintf('%s/target.txt',dir));
	W0 = load(sprintf('%s/W.%i.txt',dir,n));
	H0 = load(sprintf('%s/H.%i.txt',dir,n));
end

function create_test_data()
	target = load 'target.txt';
	W0 = load 'W.txt';
	H0 = load 'H.txt';
	save 'data.test.oct' target W0 H0;
end

% Load test data
load 'data.test.oct';

% set the save precision
save_precision(12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brunet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run Brunet algorithm
clear -all;
source('brunet/nmf.m');
function [W,H,elapsed]=run_brunet(target, W_ini, H_ini)
	rank = size(W_ini)(2)
	[tot,us,sys] = cputime();
	[W,H] = nmf(target, rank, true, W_ini, H_ini);
	[tot_e,us_e,sys_e] = cputime();
	elapsed = [tot_e,us_e,sys_e] - [tot,us,sys];
end
[W,H,elapsed] = run_brunet(target, W0, H0);
save 'ref.brunet.oct' W H elapsed;

clear -all;
source('brunet/nmf.m');
source('brunet/nmfconsensus.m');
load 'data.golub.oct';
cons = nmfconsensus(target,3,3,50,true);
cons50 = reshape(cons(3,:,:), 38,38);
save 'consensus.10.oct' cons10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNMF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear -all;
source('Benthem-Keenan-fcnnls.m');
source('Kim-Park-nmfsh_comb.m');

% Run snmf/R algorithm
function [W,H,elapsed]=run_snmfR(target, W_ini, H_ini)
	rank = size(W_ini)(2)
	[tot,us,sys] = cputime();
	[W,H,i, obj] = nmfsh_comb(target, rank, [-1 0.01], W_ini, true);
	[tot_e,us_e,sys_e] = cputime();
	elapsed = [tot_e,us_e,sys_e] - [tot,us,sys];
end
[W,H,elapsed] = run_snmfR(target, W0, H0);
save 'ref.snmfr.oct' W H elapsed;

% Run snmf/L algorithm
function [W,H,elapsed]=run_snmfL(target, W_ini, H_ini)
	rank = size(W_ini)(2)
	[tot,us,sys] = cputime();
	[Wl,Hl,i] = nmfsh_comb(target', rank, [-1 0.01], H_ini', true);
	[tot_e,us_e,sys_e] = cputime();
	elapsed = [tot_e,us_e,sys_e] - [tot,us,sys];
	W = Hl';
	H = Wl';
end
[W,H,elapsed] = run_snmfL(target, W0, H0);
save 'ref.snmfl.oct' W H elapsed;
