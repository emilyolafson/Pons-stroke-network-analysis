function [] = calculateZScores(nscans, nsess, studydir, resultsdir)
% calculateZScores: calculates the Z-score for each voxel in each patient by subtracting each patient voxel's ICC value by the mean of the controls' voxel ICC value, and diving by the standard deviation of the controls' ICC value.
% INPUT: 
%     nscans: matrix of number of scans per session (rows) per subject (columns)
%     nsess: number of sessions per subject (columns)
%     studydir: folder containing subject folders and folder for outputs.
%     resultsdir: name of folder for outputs
%
% OUTPUT: 
%     SUB*_S*_zscoreICC.mat: subject's session-specific z-scores, all GM voxels (~111000) 
%     SUB*_S*_zscoreICC.nii: subject's session-specific z-scores, in 3D (90x109x91 matrix)
%     SUBzscore: all subjects' z-scores stored in one cell array.

    load(strcat(studydir, resultsdir, 'ICC_GSR_23subjects.mat'),'ICC_GSR_cat')
    load(strcat(studydir, resultsdir, 'crtl_cat_GSR_mean.mat'), 'ctrl_cat_GSR_mean')
    load(strcat(studydir, resultsdir, 'crtl_cat_GSR_stdev.mat'), 'ctrl_cat_GSR_stdev')

    SUBzscore={} ;
    for i=1:size(nscans, 2)
        for j=1:nsess(i)        
            % load mean and stddev of controls
            SUBi_ctrlmean = ctrl_cat_GSR_mean; 
            SUBi_ctrlstdev = ctrl_cat_GSR_stdev;
        
            % subtract control mean from subject's ICC, divide by sd of controls
            %z-score calculation: (Xstrokept-mean of controls)/(standard deviation of ctls)
            SUBzscore{i}{j}=(ICC_GSR_cat{i}{j}-SUBi_ctrlmean')./SUBi_ctrlstdev';
            % save individual's z-score as a .mat
            singlesub = SUBzscore{i}{j};
            save(strcat(studydir, resultsdir, 'SUB', num2str(i), '_S', num2str(j), '_zscoreICC.mat'), 'singlesub')

            % reshape data to 3D
            GM = read_avw(strcat(studydir, 'c1referenceT1.nii'));
            GM_reshape = reshape(GM, [1 902629]);
            GM_reshape(GM_reshape > 0) = 1; %threshold GM mask
            GM_reshape(GM_reshape ==0) = 0;
            
            p=1; %counter
            zICC_GSR=[];
            zICC_GSR_3D=[];
            for z=1:size(GM_reshape,2) %starts with the gray matter mask(which has all voxels in the 3d volume.)
                % makes gray matter voxels equal to their ICC value.
                if GM_reshape(z)==1
                    zICC_GSR(z)=SUBzscore{i}{j}(p); 
                    p=p+1;
                else
                    zICC_GSR(z) =0;
                end
            end
            zICC_GSR_3D = reshape(zICC_GSR,91,109,91);
            save_avw(zICC_GSR_3D,strcat(studydir, resultsdir, 'SUB',num2str(i), '_S',num2str(j), '_zscoreICC'),'f',[2 2 2 2]);
        end
    end
    % save z-scores across subjects
    save(strcat(studydir, resultsdir, 'zICC_23subjects.mat'),'SUBzscore')
end