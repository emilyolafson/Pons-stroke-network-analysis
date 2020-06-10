clear all;

addpath([getenv('FSLDIR') '/etc/matlab']);

setenv( 'FSLDIR', '/usr/share/fsl/5.0/');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldirmpath;

% number of scans per session. row = sesssion, column = subject
numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 1 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];

% number of sessions per subject. 
nsess=[5;5;5;5;5;4;5;5;5;5;5;3;5;5;5;5;5;5;5;2;5;5;5];

studydir = '/home/emo4002/colossus_shared3/pons_sfmodelling/'
resultsdir = 'results/ICC/'

strokeptsCalculateICC(nscans, nsess, studydir, resultsdir, 'stroke_pts/')
controlsCalculateICC(ones(1, 47)+1, nsess, studydir, resultsdir, 'control_subs/');

calculateZScores(nscans, nsess,  studydir, resultsdir);

figuresdir = 'results/figures/';
disconnectivitydir = 'processing/disconnectivity/numerator_output/';

figs = [4];
k=1; %1 = cortex, 3 = cerebellum
[tst,rec,nvoxels]=makeFigures(nsess, studydir, figuresdir, resultsdir, disconnectivitydir, figs,k);
