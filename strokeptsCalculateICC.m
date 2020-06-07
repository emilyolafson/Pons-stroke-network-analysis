% Calculate ICC values om stroke patients
clear all;

addpath([getenv('FSLDIR') '/etc/matlab']);

setenv( 'FSLDIR', '/usr/share/fsl/5.0');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

% number of scans per session. row = sesssion, column = subject
numscans1_11=[4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans12_23=[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 2 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];
nscans = [numscans1_11, numscans12_23];

% number of sessions per subject. 
nsess=[5;5;5;5;5;4;5;5;5;5;5;3;5;5;5;5;5;5;5;2;5;5;5];

sub_direct='/home/emo4002/colossus_shared3/pons_sfmodelling/stroke_pts/';
%% Global signal regression and ICC calculation for concatenated scans
for i=1:23 %loop over subjects 
    %load Harvard-Oxford cortical atlas.
    % https://neurovault.org/collections/262/
    % load GM mask for global signal regression.
    GM=read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii'); 
    GM_reshape=reshape(GM, [902629 1]);
    GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
    GM_reshape(GM_reshape <=0.25) = 0;
    % get mask for regional analyses
    mask=getMask('cortex');

    %load disconnectivity files
    lesion_dc=read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/processing/get_numerator/numerator_output/SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
    lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
    
    for j=1:nsess(i) %loop over sessions

        cat_func=[]; %empty functional connectivity
        cat_outliers=[]; %empty outliers matrix
        
        cat_func = loadScans(i, j, nscans, sub_direct);
        outlierframes = scrubMotionFrames(i, j, nscans, sub_direct);
        
        
        ts=cat_func.'; %flip to <time> x <voxels>
        ts_GM = ts(:,GM_reshape==1);% <time> x <mask voxels> 
        %global signal regression    
        meants=mean(ts_GM,2); %mean across all gray matter voxels (i.e. where mask == 1)
        confounds=[meants [0; diff(meants)] outlierframes]; %confounds list
        Q=eye(size(confounds,1))-confounds*pinv(confounds);

        ts_GSR=Q*(ts_GM); %global signal regressed time series
        
        %ts_clean = high motion frames excluded.
        ts_clean_GSR=ts_GSR(sum(outlierframes,2)==0,:); %timeseries with motion frames excluded, and global signal regressed.            
        %ts_clean=ts_GM(sum(outlierframes,2)==0,:); %timeseries wiith motion frames excluded. no GSR.

        disp(strcat('saving ICC for sub:', num2str(i), ', session:', num2str(j)))
        
        %save voxelwise measure of ICC.
        ICC_GSR_cat{i}{j}=intrinsic_connectivity_contrast(ts_clean_GSR); %save timeseries into larger cell matrix
    
        % Make 3D image (for visualization)
        p=1; %counter
        for z=1:size(GM_reshape,1) %starts with the gray matter mask(which has all voxels in the 3d volume.)
            % makes gray matter voxels equal to their ICC value.
            if GM_reshape(z)==1
                ICC_GSR(z)=ICC_GSR_cat{i}{j}(p); 
                p=p+1;
            else
                ICC_GSR(z) =0;
            end
        end
        ICC_GSR_3D = reshape(ICC_GSR,91,109,91);
      %  save_avw(ICC_GSR_3D,strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/results/ICC/SUB',num2str(i), '_S',num2str(j), '_ICC'),'f',[2 2 2 2]);
    end
end

save('/home/emo4002/colossus_shared3/pons_sfmodelling/results/ICC/ICC_GSR_23subjects.mat','ICC_GSR_cat')
