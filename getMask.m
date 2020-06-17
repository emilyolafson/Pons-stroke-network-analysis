function mask = getMask(region, studydir)
% OPTIONS
% region: cortex, cerebellum, left primary motor cortex (left m1), right primary motor cortex (m2), none
    if (strcmp(region,'left m1') | strcmp(region, 'right m1'))
        cort_atlas = read_avw(strcat('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-2mm.nii.gz')); %copied from /conn/utils/surf/c1referenceT1.nii. 91x109x91
        cort_struct = parseXML(strcat('/usr/share/fsl/5.0/data/atlases/HarvardOxford-Cortical.xml'));
        subcort_atlas = read_avw(strcat('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr0-2mm.nii.gz')); %copied from /conn/utils/surf/c1referenceT1.nii. 91x109x91
        subcort = reshape(subcort_atlas,[],1); %flattened 1D <voxel> matrix with probability [0,1] being gray matter
        cort = reshape(cort_atlas,[],1); %flattened 1D <voxel> matrix with probability [0,1] being gray matter

        template=zeros(1,902629);
        if strcmp(region, 'left m1')
            template(cort== 7) = 1;
            template(subcort ~= 2) = 0;
             mask = reshape(template, [91 109 91]);
        end
        if strcmp(region, 'right m1')
            template(cort_atlas == 7) = 1;
            template(subcort_atlas ~= 13)= 0;
            mask = reshape(template, [91 109 91]);

        end
        
    end
    
    if (strcmp(region, 'cortex') | strcmp(region,'brainstem'))
        subcort_atlas = read_avw(strcat('/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz')); %copied from /conn/utils/surf/c1referenceT1.nii. 91x109x91
        subcort_struct = parseXML(strcat('/usr/share/fsl/5.0/data/atlases/HarvardOxford-Subcortical.xml'));
        idx = 1;
        for i = 2:2:42
            subcort_region{idx} = subcort_struct.Children(4).Children(i).Children.Data;
            idx = idx+1;
        end
        subcort = reshape(subcort_atlas,[],1); %flattened 1D <voxel> matrix with probability [0,1] being gray matter
        template = zeros(1, 902629);
        if strcmp(region, 'cortex')
            template(subcort == 13 | subcort == 2) = 1;
            mask = reshape(template, [91 109 91]);
        end
        if strcmp(region, 'brainstem')
            template(subcort == 8) = 1;
            mask = reshape(template, [91 109 91]);
        end
    end

    if (strcmp(region,'cerebellum') | strcmp(region,'r cerebellum') | strcmp(region,'l cerebellum')) 
        cere_atlas = read_avw(strcat('/usr/share/fsl/5.0/data/atlases/Cerebellum/Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz')); %copied from /conn/utils/surf/c1referenceT1.nii. 91x109x91
        cere_struct = parseXML(strcat('/usr/share/fsl/5.0/data/atlases/Cerebellum_MNIfnirt.xml'));
        idx = 1;
        for i = 2:2:56
            cere_region{idx} = cere_struct.Children(4).Children(i).Children.Data;
            idx = idx+1;
        end
        cere = reshape(cere_atlas,[],1); %flattened 1D <voxel> matrix with probability [0,1] being gray matter
        template = zeros(1, 902629);
        template(cere > 0) = 1;
        mask = reshape(template, [91 109 91]);
    end
    if (strcmp(region, 'none'))
        GM = read_avw(strcat(studydir, 'c1referenceT1.nii')); 
        GM_reshape = reshape(GM, [902629 1]);
        mask= GM_reshape;
    end

end
