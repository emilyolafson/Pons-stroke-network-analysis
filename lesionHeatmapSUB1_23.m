% calculate lesion mask heatmap
clear all;
if(isempty(which('save_avw'))) %add FSL functions to matlab path.
        addpath(['/usr/share/fsl/5.0/' '/etc/matlab']);
end

%% load lesionMask data & create sum
%cell array of lesionMask data.
for i =1:23
    subject=strcat('SUB', num2str(i))
    lesionmask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/stroke_pts/',subject,'/lesionMask/wc_',subject,'_lesion.nii.gz')); %load MNI 2x2x2mm mask data
    lesionmask{i}=reshape(lesionmask{i},[],size(lesionmask{i},4)); %flattened 2D matrix that is <voxels>
end

%create <subject x voxels matrix>
all_lesions = [lesionmask{1}'; lesionmask{2}'; lesionmask{3}'; lesionmask{4}'; lesionmask{5}'; lesionmask{6}'; lesionmask{7}'; lesionmask{8}'; lesionmask{9}'; lesionmask{10}';lesionmask{11}'; lesionmask{12}'; lesionmask{13}'; lesionmask{14}'; lesionmask{15}'; lesionmask{16}'; lesionmask{17}'; lesionmask{18}'; lesionmask{19}'; lesionmask{20}';lesionmask{21}'; lesionmask{22}'; lesionmask{23}'];

% take sum across subjects
lesion_sum = sum(all_lesions);

% reshape to 3D voxel volume
reshape_sum = reshape(lesion_sum, 91,109,91);

%
maxoverlap=max(max(max(reshape_sum)));

% save as .nii file to display
save_avw(reshape_sum,strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/results/subjectSum_lesionMask'),'f',[2 2 2 2]);

