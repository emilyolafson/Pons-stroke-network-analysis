%Calcualte ICC across control subjects & save mean & std dev for each
%voxel.
% change test. 
clear all;

addpath([getenv('FSLDIR') '/etc/matlab']);
setenv( 'FSLDIR', '/usr/share/fsl/5.0');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
clear fsldir fsldirmpath;

sub_direct='/home/emo4002/colossus_shared3/pons_sfmodelling/control_subs/';
%% Global signal regression and ICC calculation for concatenated scans
for i=24:47 %loop over subjects   
    %load Harvard-Oxford cortical atlas.
    % https://neurovault.org/collections/262/
    
     % load GM mask for global signal regression.
    GM=read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii'); 
    GM_reshape=reshape(GM, [902629 1]);
    GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
    GM_reshape(GM_reshape <=0.25) = 0;
    
    % get mask for regional analyses
    mask=getMask('cortex');
    
    j=1;
    
    cat_func=[]; %empty functional connectivity
    cat_outliers=[]; %empty outliers matrix
    nscans=ones(1, 24)+1;
    cat_func = loadScans(i, j, nscans, sub_direct);
    outlierframes = scrubMotionFrames(i, j, nscans);

    ts=cat_func.'; %flip to <time> x <voxels>
    ts_GM = ts(:,GM==1);% <time> x <mask voxels> 

    %global signal regression    
    meants=mean(ts_GM,2); %mean across all gray matter voxels (i.e. where mask == 1)
    confounds=[meants [0; diff(meants)] outlierframes]; %confounds list
    Q=eye(size(confounds,1))-confounds*pinv(confounds);
    ts_GSR=Q*(ts_GM); %global signal regressed time series

    %ts_clean = high motion frames excluded.
    ts_clean_GSR=ts_GSR(sum(outlierframes,2)==0,:); %timeseries with motion frames excluded, and global signal regressed.            
    %ts_clean=ts_GM(sum(outlierframes,2)==0,:); %timeseries wiith motion frames excluded. no GSR.

    disp(strcat('saving ICC for ', subject, ' session ', num2str(j)))

    %save voxelwise measure of ICC.
    ICC_GSR_cat_controls{i}=intrinsic_connectivity_contrast(ts_clean_GSR); %save timeseries into larger cell matrix

    % Make 3D image (for visualization)
    p=1; %counter
    for z=1:size(GM,1) %starts with the gray matter mask(which has all voxels in the 3d volume.)
        % makes gray matter voxels equal to their ICC value.
        if GM(z)==1
            ICC_GSR(z)=ICC_GSR_cat_controls{i}(p); 
            p=p+1;
        else
            ICC_GSR(z) =0;
        end
    end
    ICC_GSR_3D = reshape(ICC_GSR,91,109,91);
    save_avw(ICC_GSR_3D,strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/results/ICC/',subject, '_S',num2str(j), '_ICC'),'f',[2 2 2 2]);
end

save('/home/emo4002/colossus_shared3/pons_sfmodelling/results/ICC/ICC_GSR_controlsubjects.mat','ICC_GSR_cat_controls')

%% Calculate mean and standard deviation of ICC across control subjects - from concatenated scans.

%load('/home/emo4002/colossus_shared3/pons_sfmodelling/results/ICC/ICC_GSR_controlsubjects.mat', 'ICC_GSR_cat_controls');
cat=[];
cat_GSR=[];
for i=24:47
    subject = strcat('SUB', num2str(i));
    GM=read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii'); 
    GM_reshape=reshape(GM, [902629 1]);
    GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
    GM_reshape(GM_reshape <=0.25) = 0;
    
    ICC_GSR_cat_controls = reshape(ICC_GSR_final_cat{i}, [], size(ICC_GSR_final_cat{i}, 4));

    cat=[cat; ICC_cat'];
    cat_GSR=[cat_GSR;ICC_GSR_cat_controls']; %the img is stored in the 2nd cell. add each new subj as a row.
end

ctrl_cat_GSR_mean = mean(cat_GSR, 'omitnan');
ctrl_cat_GSR_stdev= std(cat_GSR,'omitnan');
ctrl_cat_mean= mean(cat, 'omitnan');
ctrl_cat_stdev= std(cat,'omitnan');

save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_GSR_mean.mat', 'ctrl_cat_GSR_mean')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_mean.mat', 'ctrl_cat_mean')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_GSR_stdev.mat', 'ctrl_cat_GSR_stdev')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_stdev.mat', 'ctrl_cat_stdev')
