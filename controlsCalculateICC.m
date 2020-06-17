function [] = controlsCalculateICC(nscans, nsess, studydir, resultsdir, controldir)
    % controlsCalculateICC: Calculate ICC values in 24 control patients
    % INPUT: 
    %     nscans: matrix of number of scans per session (rows) per subject (columns)
    %     nsess: number of sessions per subject (columns)
    %     studydir: folder containing subject folders and folder for outputs.
    %     resultsdir: name of folder for outputs
    %     controldir: name of control subjects folder.
    %
    % OUTPUT: 
    %     SUB*_S*_ICC.nii.gz: 3D rendering of ICC values in 2mm MNI space (91x109x91 voxels) 
    %     ICC_GSR_cat.mat: the global signal regressed ICC values across all voxels in the gray matter in a cell array, where each cell contains the ICC values for each session saved separately.

    

    addpath([getenv('FSLDIR') '/etc/matlab']);
    setenv( 'FSLDIR', '/usr/share/fsl/5.0');
    fsldir = getenv('FSLDIR');
    fsldirmpath = sprintf('%s/etc/matlab',fsldir);
    path(path, fsldirmpath);
    clear fsldir fsldirmpath;

    %% Global signal regression and ICC calculation for concatenated scans
    for i=24:47 %loop over controls subjects   
    
        % load GM mask for global signal regression.
        GM=read_avw(strcat(studydir, '/c1referenceT1.nii')); 
        GM_reshape=reshape(GM, [902629 1]);
        GM_reshape(GM_reshape > 0) = 1; %threshold GM mask
        GM_reshape(GM_reshape ==0) = 0;
      
        cat_func=[]; %empty functional connectivity
        cat_outliers=[]; %empty outliers matrix
    
        j=1;
    
        sub_direct = strcat(studydir, controldir);
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

        disp(strcat('saving ICC for ', num2str(i), ' session ', num2str(j)))

        %save voxelwise measure of ICC.
        ICC_GSR_cat_controls{i-23}=intrinsic_connectivity_contrast(ts_clean_GSR); %save timeseries into larger cell matrix

        % Make 3D image (for visualization)
        p=1; %counter
        ICC_GSR=[];
        ICC_GSR_3D=[];
        for z=1:size(GM_reshape,1) %starts with the gray matter mask(which has all voxels in the 3d volume.)
            % makes gray matter voxels equal to their ICC value.
            if GM_reshape(z)==1
                ICC_GSR(z)=ICC_GSR_cat_controls{i-23}(p); 
                p=p+1;
            else
                ICC_GSR(z) =0;
            end
        end
        ICC_GSR_3D = reshape(ICC_GSR,91,109,91);
        save_avw(ICC_GSR_3D,strcat(studydir, resultsdir, 'SUB',num2str(i), '_ICC'),'f',[2 2 2 2]);
    end

    save(strcat(studydir, resultsdir, 'ICC_GSR_controlsubjects.mat'),'ICC_GSR_cat_controls')

    %% Calculate mean and standard deviation of ICC across control subjects - from concatenated scans.

    %load(studydir, resultsdir, 'ICC_GSR_controlsubjects.mat', 'ICC_GSR_cat_controls');
    cat_GSR=[];
    for i=24:47
        cat_GSR = [cat_GSR, ICC_GSR_cat_controls{i-23}'];
    end
    ctrl_cat_GSR_mean = mean(cat_GSR, 2);
    ctrl_cat_GSR_stdev= std(cat_GSR')';

    save(strcat(studydir, resultsdir, 'crtl_cat_GSR_mean.mat'), 'ctrl_cat_GSR_mean')
    save(strcat(studydir, resultsdir, 'crtl_cat_GSR_stdev.mat'), 'ctrl_cat_GSR_stdev')
end