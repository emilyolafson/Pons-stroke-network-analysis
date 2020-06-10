%average lesionTract images
clear all;
if(isempty(which('save_avw'))) %add FSL functions to matlab path.
        addpath([getenv('FSLDIR') '/etc/matlab']);
end

%cell array of lesionTract data.
for i =1:11
    subject=strcat('SUB', num2str(i))
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz')); %load lesiontract data
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>

    %mask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionMask/wc_',subject,'_lesion.nii.gz')); %load lesion Mask data
    %mask{i}=reshape(mask{i},[],size(mask{i},4)); %flattened 2D matrix that is <voxels>

end

all_lesions = [lesion_dc{1}'; lesion_dc{2}'; lesion_dc{3}'; lesion_dc{4}'; lesion_dc{5}'; lesion_dc{7}'; lesion_dc{7}'; lesion_dc{8}'; lesion_dc{9}'; lesion_dc{10}'];
%all_masks = [mask{1}'; mask{2}'; mask{3}'; mask{4}'; mask{5}'; mask{7}'; mask{7}'; mask{8}'; mask{9}'; mask{10}'];

%average values across subjects.
lesionAvg = mean(all_lesions);
%maskAvg = sum(all_masks);

reshape_avg = reshape(lesionAvg, 91,109,91);
%reshape_avg_mask = reshape(maskAvg, 91,109,91);

save_avw(reshape_avg,strcat('/home/emo4002/colossus_shared3/stroke_pts/SavgLesionTract'),'f',[2 2 2 2]);
%save_avw(reshape_avg_mask,strcat('/home/emo4002/colossus_shared3/stroke_pts/SavgMask'),'f',[2 2 2 2]);
