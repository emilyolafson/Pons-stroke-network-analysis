clear all;

addpath([getenv('FSLDIR') '/etc/matlab']);

setenv( 'FSLDIR', '/usr/share/fsl/5.0');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldirmpath;

% number of scans per session. row = sesssion, column = subject
numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 2 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];

% number of sessions per subject. 
nsess=[5;5;5;5;5;4;5;5;5;5;5;3;5;5;5;5;5;5;5;2;5;5;5];

studydir = '/home/emo4002/colossus_shared3/pons_sfmodelling/'
resultsdir = 'results/ICC/'
strokedir='stroke_pts/'
controldir='control_subs/'
strokeptsCalculateICC(nscans, nsess, studydir, resultsdir, strokedir)
controlsCalculateICC(ones(1, 47)+1, nsess, studydir, resultsdir, 'control_subs/');

calculateZScores(nscans, nsess,  studydir, resultsdir);

figuresdir = 'results/figures/';
disconnectivitydir = 'processing/disconnectivity/numerator_output/';

figs = [2,8];

k=1; %1 = cortex, 2 = cerebellum, 3 = left motor cortex, 4 = right motor cortex

[tst,rec,nvoxels]=makeFigures(nsess, studydir, strokedir, figuresdir, resultsdir, disconnectivitydir, figs,k);
    % Figure 1 - Boxplots of z-scores of ICC of voxels connected vs unconnected to lesion
    % Figure 2 - Change in motor scores correlated with the change in z-score ICC from the last scan to the first scan 
    % Figure 3 - Correlation betwen inter-session change in ICC & inter-session change in Fugl-Meyer (i.e., Session 2 vs Session 1, Session 3 vs Session 2)
    % Figure 4 - Change in ICC in areas structurally connected to the lesion for each subject (histogram of values @ first & last session)
    % Figure 5 - Same as Figure 3, but each inter-session change is plotted separately.
    % Figure 6 - plotting motor recovery information (Fugl Meyer)
    % Figure 7 - plotting motor recovery information (Grip Strength)
    % Figure 8 - Comparison between NeMo2 voxelwise scores and LEAD-DBS scores
    % Figure 9 - Voxelwise correlation between z-score ICC and structural disconnectivity to lesion (NeMo) 
    % Figure 10 - Correlation between ICC scores in affected hemisphere and grip strength in affected hand
    % Figure 11 - Correlation between average LEAD-DBS output and average NeMo2 oiutput (across 11 subjects)
