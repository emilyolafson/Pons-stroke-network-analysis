%Remove high-motion frames, perform GSR (or not), calculate ICC.
clear all;
lesion_dc={{},{},{},{},{}}
graymask={{},{},{},{},{}}
sessions ={{},{},{},{},{}}
scans ={{},{},{},{},{}}
sessions={scans,scans,scans,scans,scans}
func ={sessions,sessions,sessions,sessions,sessions,sessions,sessions,sessions,sessions,sessions,sessions}
ICC_GSR=func;
ICC=func;
handle_gm_GSR=func;
handle_gm=func;
Cupper_gm=func;
Cupper_gsr= func;

if(isempty(which('save_avw')))
    addpath([getenv('FSLDIR') '/etc/matlab']);
end

numscans = [4 3 3 3 3 3 3 3 3 3 3;4 3 3 3 3 3 3 4 3 3 3;3 3 3 3 3 3 3 3 3 3 3;4 4 3 3 3 3 3 3 3 3 3 ; 3 3 3 3 3 0 3 3 3 3 3];
numscans2 =[3 3 3 3 3 2 2 2 2 2 2 2; 3 3 3 4 3 2 2 2 2 2 2 2 ;3 3 3 3 3 2 2 1 0 2 2 2;0 3 3 3 3 2 2 2 0 2 2 2;0 3 4 3 3 2 2 2 0 2 2 2];

nsess=[5;5;5;5;5;4;5;5;5;5;5];
%% Global signal regression and ICC calculation.
% apply GM mask to lesionTract data (multiply binary gray matter mask)
% apply GM mask to EPI data (take subset where 
    % this is because the nd/or ICC function will computer the R^2 correlation
    % between a given voxel and all other voxels. Thus having many voxels
    % with value 0 will artificially drop ICC values. 
% perform global signal regression (GSR)
% generate plots of correlation  coefficients before and after GSR.
% compute ICC values on gray matter voxels
% generate plots of ICC values with and without GSR.
for i =19:44
    if i==36
        continue
    end
    if i==39
        continue
    end
    
    subject = strcat('SUB', num2str(i));
    graymask{i}=readavw(strcat('/home/emo4002/colossus_shared3/c1referenceT1.nii')
   % graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5

   % lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
   % lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
   % lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels>
    for j = 1:nsess(i)
        for k =1:2
            disp(strcat('loading functional scans for ', subject, ' session ', num2str(j), ' scan ', num2str(k)))
            func{i}{j}{k}=read_avw(strcat('/home/emo4002/colossus_shared3/control_subs/',subject,'/func/denoise_swaufunc',num2str(k),'.nii'));
            func{i}{j}{k}=reshape(func{i}{j}{k},[],size(func{i}{j}{k},4)); %flattened 2D matrix that is <voxels>
            %func{i}{j}{k}=func{i}{j}{k} %same but now <voxels> x <time> where voxels not in the mask are set to zero
            ts=func{i}{j}{k}.'; %flip to <time> x <voxels>
            outliers=load(strcat('/home/emo4002/colossus_shared3/control_subs/',subject,'/func/art_regression_outliers_aufunc_',num2str(k),'.mat'));
            outlierframes = outliers.mat;
            outlierframes = outlierframes(:, sum(outlierframes,1) > 0); %
                 
            motionframes =find(sum(outlierframes,2));
            motionframes=unique([motionframes; motionframes+1]); %add frame after high-motion frame (frame before is already flagged as an outlier by CONN).
            motionframes = motionframes(motionframes<=124); %cut off frame 125 if it gets added by step above.
            outlierframes=zeros(124, size(motionframes,1)); %outlierframes = empty matrix to be populated below.
            %set cells to 1 if they are outliers (one outlier per row/column) (adds frame after ART-identified frame - the one
            %before has already been set to 1, by ART)
            for l=1:size(motionframes,1)
                outlierframes(motionframes(l),l) = 1;
            end
            disp(strcat('confound regression for ', subject, ' session ', num2str(j), ' scan ', num2str(k)))
            
            % confound regression
            ts_GM = ts(:,graymask{i}==1);% time x gray matter voxels 
            meants=mean(ts_GM,2); %mean across all gray matter voxels (i.e. where mask == 1)
            confounds=[meants [0; diff(meants)] outlierframes]; %confounds list
            Q=eye(size(confounds,1))-confounds*pinv(confounds);
            ts_GSR=Q*(ts_GM); %global signal regressed time series
            
            [~,ridx]=sort(rand(size(ts_GM,2),1)); 
            rand1000=ridx(1:1000); %sample and store 1000 random voxels

           % figure(1) %heatmap of voxels without GSR
            handle_gm{i}{j}{k}=heatmap(ts_GM(sum(outlierframes,2)==0,rand1000)');

           % figure(2) %heatmap with GSR
            handle_gm_GSR{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0,rand1000)');

              
            %correlation coefficients of the 1000 voxels 
            %in non-GSR gray matter
            C_gm=corr(ts_GM(:,rand1000)); %should be 1000x1000
            trimask=triu(ones(size(C_gm)),1)>0;
            Cupper_gm{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
            % in GSR gray matter.
            C_gsr=corr(ts_GSR(:,rand1000)); %should be 1000x1000
            trimask=triu(ones(size(C_gsr)),1)>0;
            Cupper_gsr{i}{j}{k}=C_gsr(trimask); %should be %1000*1000*/2-1000
            
          %  figure(4) %difference in voxel values before and after gm
           % handle_gm_GSR_DIFF{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0, rand1000)'-ts_GM(sum(outlierframes,2)==0,rand1000)');
            
            ts_clean_GSR=ts_GSR(sum(outlierframes,2)==0,:); %timeseries with motion frames excluded, and global signal regressed.            
            ts_clean=ts_GM(sum(outlierframes,2)==0,:); %timeseries wiith motion frames excluded. no GSR.

            disp(strcat('saving ICC for ', subject, ' session ', num2str(j), ' scan ', num2str(k)))

            ICC_GSR{i}{j}{k} =intrinsic_connectivity_contrast(ts_clean_GSR); %save timeseries into larger cell matrix
            ICC{i}{j}{k} = intrinsic_connectivity_contrast(ts_clean); % "
 
            % MAKE FULL 3D IMAGE (for displaying)
            p=1; %counter
            a = graymask{i};
            for z=1:size(graymask{i},1) %starts with the gray matter mask(which has all voxels in the 3d volume.)
                % makes gray matter voxels equal to their ICC value.
                if a(z)~=0
                    full_ICC(z)=ICC{i}{j}{k}(p); %first GM value of the volume is 
                    full_ICC_GSR(z)=ICC_GSR{i}{j}{k}(p); 
                    p=p+1;
                else
                    full_ICC(z) = 0;
                    full_ICC_GSR(z) =0;
                end
            end
        %   ICC_final{i}{j}{k} = reshape(full_ICC,91,109,91);
            ICC_GSR_final{i}{j}{k} = reshape(full_ICC_GSR,91,109,91);
            %save_avw(test2,'/home/emo4002/colossus_shared3/pons_sfmodelling/20test','f',[2 2 2 2]);
            save_avw(ICC_GSR_final{i}{j}{k},strcat('/home/emo4002/colossus_shared3/control_subs/',subject, '/func/S',num2str(j), '/ICC_GSR',subject, num2str(k)),'f',[2 2 2 2]);
          % save_avw(ICC_final{i}{j}{k},strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject, '/func/S',num2str(j), '/ICC',subject, num2str(k)),'f',[2 2 2 2]);
        end
    end
end

save('/home/emo4002/colossus_shared3/control_subs/ICC_GSR_controls.mat','ICC_GSR')
%save('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_controls.mat','ICC')

%% Global signal regression and ICC calculation - CONCATENATED SCANS
for i =19:44
    if i==36
        continue
    end
    if i==39
        continue
    end
    subject = strcat('SUB', num2str(i));
    graymask{i}=read_avw(strcat('/home/emo4002/colossus_shared3/c1referenceT1.nii'))
    % graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel

    %lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
   % lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
   % lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels>
    for j = 1
        cat_func=[];
        cat_outliers=[];
        for k =1:2
            disp(strcat('loading functional scans for ', subject, ' session ', num2str(j), ' scan ', num2str(k)))
            disp(num2str(i))
            func{i}{j}{k}=read_avw(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/func/denoise_swaufunc',num2str(k),'.nii'));
            func{i}{j}{k}=reshape(func{i}{j}{k},[],size(func{i}{j}{k},4)); %flattened 2D matrix that is <voxels>
            %func{i}{j}{k}=func{i}{j}{k} %same but now <voxels> x <time> where voxels not in the mask are set to zero
            cat_func = horzcat(cat_func, func{i}{j}{k}); % concatenate scans horizontally.
        end
        outliers=load(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/func/art_regression_outliers_aufunc1.mat'));
        outlierframes = outliers.R;
        outlierframes = outlierframes(:, sum(outlierframes,1) > 0); %


        motionframes =find(sum(outlierframes,2));
        motionframes=unique([motionframes; motionframes+1]); %add frame after high-motion frame (frame before is already flagged as an outlier by CONN).
        motionframes = motionframes(motionframes<=size(outlierframes,1)); %cut off frame 125 if added.
        outlierframes=zeros(size(outlierframes,1), size(motionframes,1)); %outlierframes = empty matrix to be populated below.
        %set cells to 1 if they are outliers (one outlier per
        %row/column) (adds frame after ART-identified frame - the one
        %before has already been set to 1, by ART)

        for l=1:size(motionframes,1)
            outlierframes(motionframes(l),l) = 1;
        end

        ts=cat_func.'; %flip to <time> x <voxels>
   
            
        disp(strcat('confound regression for ', subject, ' session ', num2str(j)))

        % confound regression
        ts_GM = ts(:,graymask{i}>0.5);% time x gray matter voxels 

        meants=mean(ts_GM,2); %mean across all gray matter voxels (i.e. where mask == 1)
        confounds=[meants [0; diff(meants)] outlierframes]; %confounds list
        Q=eye(size(confounds,1))-confounds*pinv(confounds);
        ts_GSR=Q*(ts_GM); %global signal regressed time series

        [~,ridx]=sort(rand(size(ts_GM,2),1)); 
        rand1000=ridx(1:1000); %sample and store 1000 random voxels

       % figure(1) %heatmap of voxels without GSR
       % handle_gm_cat{i}{j}{k}=heatmap(ts_GM(sum(outlierframes,2)==0,rand1000)');

       % figure(2) %heatmap with GSR
        handle_gm_GSR_cat{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0,rand1000)');

        %correlation coefficients of the 1000 voxels 
        %in non-GSR gray matter
      %  C_gm=corr(ts_GM(:,rand1000)); %should be 1000x1000
      %  trimask=triu(ones(size(C_gm)),1)>0;
       % Cupper_gm_cat{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
        % in GSR gray matter.
        C_gsr=corr(ts_GSR(:,rand1000)); %should be 1000x1000
        trimask=triu(ones(size(C_gsr)),1)>0;
        Cupper_gsr_cat{i}{j}{k}=C_gsr(trimask); %should be %1000*1000*/2-1000

        handle_gm_GSR_DIFF_cat{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0, rand1000)'-ts_GM(sum(outlierframes,2)==0,rand1000)');

        ts_clean_GSR=ts_GSR(sum(outlierframes,2)==0,:); %timeseries with motion frames excluded, and global signal regressed.            
      %  ts_clean=ts_GM(sum(outlierframes,2)==0,:); %timeseries wiith motion frames excluded. no GSR.

        disp(strcat('saving ICC for ', subject, ' session ', num2str(j)))

        ICC_GSR_controls_cat{i}{j}{k} =intrinsic_connectivity_contrast(ts_clean_GSR); %save timeseries into larger cell matrix
      %  ICC_controls_cat{i}{j}{k} = intrinsic_connectivity_contrast(ts_clean); % "
        %   figure(5)
        %   ksdensity(ICC_GSR{i}{j}{k})
        %   title('ICC values after GSR')
        %   figure(6)
        %   ksdensity(ICC{i}{j}{k})
        %   title('ICC values before GSR')

        % MAKE FULL 3D IMAGE
        p=1; %counter
        a = graymask{i};
        %for z=1:size(graymask{i},1) %starts with the gray matter mask(which has all voxels in the 3d volume.)
            % makes gray matter voxels equal to their ICC value.
      %      if a(z)~=0
     %     %      full_ICC(z)=ICC_controls_cat{i}{j}{k}(p); %first GM value of the volume is 
                full_ICC_GSR(z)=ICC_GSR_controls_cat{i}{j}{k}(p); 
     %           p=p+1;
      %      else
              %  full_ICC(z) = 0;
     %           full_ICC_GSR(z) =0;
     %       end
     %   end
      %  ICC_final_cat{i}{j}{k} = reshape(full_ICC,91,109,91);
        ICC_GSR_final_cat{i}{j}{k} = reshape(full_ICC_GSR,91,109,91);
        %save_avw(test2,'/home/emo4002/colossus_shared3/pons_sfmodelling/20test','f',[2 2 2 2]);
        save_avw(ICC_GSR_final_cat{i}{j}{k},strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject, '/func/ICC_GSR_cat', subject, 'S', num2str(j)),'f',[2 2 2 2]);
      %  save_avw(ICC_final_cat{i}{j}{k},strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject, '/func/ICC_cat', subject, 'S', num2str(j)),'f',[2 2 2 2]);
    end
end

%save('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_final_cat_controls.mat', 'ICC_final_cat')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_GSR_final_cat_controls.mat', 'ICC_GSR_final_cat')

save('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_GSR_controls_cat.mat','ICC_GSR_controls_cat')
%save('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_controls_cat.mat','ICC_controls_cat')

%% Calculate mean and standard deviation of ICC across control subjects - from concatenated scans.

load('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_final_cat_controls.mat', 'ICC_final_cat');
load('/home/emo4002/colossus_shared3/pons_sfmodelling/ICC_GSR_final_cat_controls.mat', 'ICC_GSR_final_cat');
cat=[];
cat_GSR=[];
for i=19:44
    
    if i==36
        continue
    end
    if i==39
        continue
    end
    subject = strcat('SUB', num2str(i));
    graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/CTL_prep/cleanData/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel

    disp(num2str(i))
    ICC_cat = reshape(ICC_final_cat{i}{1}{2}, [], size(ICC_final_cat{i}{1}{2}, 4)); %load so that you have <voxels>
    ICC_GSR_cat = reshape(ICC_GSR_final_cat{i}{1}{2}, [], size(ICC_GSR_final_cat{i}{1}{2}, 4));
    
    ICC_cat(ICC_cat(graymask{i}==1)==0)=NaN; %set values inside the gray mask that have an ICC of zero to be "NaN". This ideally eliminates any problems with the masks not fully lining up between subjects
    ICC_GSR_cat(ICC_GSR_cat(graymask{i}==1)==0)=NaN;

    cat=[cat; ICC_cat'];
    cat_GSR=[cat_GSR;ICC_GSR_cat']; %the img is stored in the 2nd cell. add each new subj as a row.
end

ctrl_cat_GSR_mean = mean(cat_GSR, 'omitnan');
ctrl_cat_GSR_stdev= std(cat_GSR,'omitnan');
ctrl_cat_mean= mean(cat, 'omitnan');
ctrl_cat_stdev= std(cat,'omitnan');

save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_GSR_mean.mat', 'ctrl_cat_GSR_mean')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_mean.mat', 'ctrl_cat_mean')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_GSR_stdev.mat', 'ctrl_cat_GSR_stdev')
save('/home/emo4002/colossus_shared3/pons_sfmodelling/crtl_cat_stdev.mat', 'ctrl_cat_stdev')

%% correlation between ICC values before and after GSR with averaging vs concatenation
p=1; %subplot counter

%store correlation coefficients
    %rows: subjects
    %columns: sessions (longitudinal time points)

%ICC values concatenating vs averaging
conc_vs_avg_GSR=[];
conc_vs_avg_noGSR=[];

%ICC values before vs after GSR
before_vs_after_conc=[]; 
before_vs_after_avg=[];

%loop over subjects, store correlations between ICC conc/avg and GSR
%before/after.
for i=1:5
    p=1;
    subject = strcat('SUB', num2str(i));
    graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel

    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels>
    for j=1:nsess(i)
        scans_gsr=[];
        scans=[];
        %compute average ICC per session 
        for k=1:numscans(j,i)
            scans_gsr=[scans_gsr; reshape(ICC_GSR_final{i}{j}{k},[],1)']; %masks out things not in the gray matter for plotting.
            scans =[scans; reshape(ICC_final{i}{j}{k},[],1)'];
        end
        %gsr ICC values
        icc_avg_GSR = mean(scans_gsr);
        icc_concat_GSR = reshape(ICC_GSR_final_cat{i}{j}{k},[],1);

        %non gsr ICC values
        icc_avg = mean(scans);
        icc_concat = reshape(ICC_final_cat{i}{j}{k}, [],1);
        
        a=corrcoef(icc_concat, icc_avg); %cconcatenated vs average - no gsr
        b=corrcoef(icc_concat_GSR, icc_avg_GSR); % concatenated vs average - gsr
        c=corrcoef(icc_avg, icc_avg_GSR); % before vs after gsr - average
        d=corrcoef(icc_concat, icc_concat_GSR); %before vs after gsr - concatenated
        
        conc_vs_avg_noGSR(i,j)= a(1,2);
        conc_vs_avg_GSR(i,j)=b(1,2);
        before_vs_after_avg(i,j)=c(1,2);
        before_vs_after_conc(i,j)=d(1,2);
        
        lesion = lesion_dc{i}(graymask{i}==1);
        figure(i)
        subplot(1,nsess(i),p);
        %scatterplot: ICC for each voxel from averaging (yaxis) vs from
        %concatenating (yaxis)
        %no global signal regression
        %scatter(icc_concat, icc_avg, '.r', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeAlpha', 0.1);
        scatter(icc_avg, icc_avg_GSR,'.g', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeAlpha', 0.1)
        hold on
        %global signal regressed
        %scatter(icc_concat_GSR, icc_avg_GSR, '.b', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeAlpha', 0.1);
        scatter(icc_concat, icc_concat_GSR, '.b', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeAlpha', 0.1)
        hold on;
        refline(1)
        xlabel('ICC before GSR')
        ylabel('ICC after GSR')
        xlim([0 0.1])
        ylim([0 0.1])
        legend(['Averaging'],['Concatenating'])
        title(['Subject ', num2str(i), 'Session ', num2str(j)])
        p=p+1;
        set(gcf, 'Position', [10 10 2000 300]);
     
    end
end

figure(2)
subplot(1,2,2)
boxplot(before_vs_after_conc);
ylim([0.95 1])
xticklabels({'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'})
title('Correlation between ICC values before and after GSR - Concatenating scans')
ylabel('Correlation Coefficients')
xlabel('Longitundinal time point')

subplot(1,2,1)
boxplot(before_vs_after_avg);
ylim([0.95 1])
xticklabels({'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'})
ylabel('Correlation Coefficients')
title('Correlation between ICC values before and after GSR - Averaging scans')
xlabel('Longitundinal time point')

hold on;
set(gcf, 'Position', [10 10 2000 1000])

%% within-voxel effect on ICC of averaging vs concatenating scans
figure(3)
subplot(1,2,2)
boxplot(conc_vs_avg_GSR);
ylim([0.95 1])
xticklabels({'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'})
title('Correlation between ICC values averaging vs concatenating scans - with GSR')
ylabel('Correlation Coefficients')
xlabel('Longitundinal time point')

subplot(1,2,1)
boxplot(conc_vs_avg_noGSR);
ylim([0.95 1])
xticklabels({'Session 1', 'Session 2', 'Session 3', 'Session 4', 'Session 5'})
ylabel('Correlation Coefficients')
title('Correlation between ICC values averaging vs concatenating scans - no GSR')
xlabel('Longitundinal time point')

hold on;
set(gcf, 'Position', [10 10 2000 1000])


%% compute average ICC across all subjects at time point 1. (not concatenated)
all_ICC={};
all_ICC_GSR={};
diff_ICC= [];
sess_avg=sessions;

for i = 1:11
    for j=1:nsess(i)
         n=numscans(j,i)
         all_ICC(1:n)={{}};
         all_ICC_GSR(1:n)={{}};
        for k=1:n
            all_ICC{k} = [all_ICC{k} ; reshape(ICC_final{i}{j}{k},[],1)'];
            all_ICC_GSR{k} = [all_ICC_GSR{k} ; reshape(ICC_GSR_final{i}{j}{k},[],1)'];
        end
        All = [all_ICC{:}];
        All_GSR = [all_ICC_GSR{:}];

        if n==3
           all = [All{n};All{n-1};All{n-2}];
           all_gsr = [All_GSR{n};All_GSR{n-1};All_GSR{n-2}];
        end
        
        if n==4
           all = [All{n};All{n-1};All{n-2};All{n-3}];
           all_gsr = [All_GSR{n};All_GSR{n-1};All_GSR{n-2};All_GSR{n-3}];

        end     
        sess_avg{i}{j}=mean(all,1);
        sess_avg_gsr{i}{j}=mean(all_gsr,1);

    end
    
    if i==6
         sub_avg{i}=mean([sess_avg{i}{1};sess_avg{i}{2};sess_avg{i}{3};sess_avg{i}{4}])
         sub_avg_gsr{i}=mean([sess_avg_gsr{i}{1};sess_avg_gsr{i}{2};sess_avg_gsr{i}{3};sess_avg_gsr{i}{4}])

    else 
         sub_avg{i}=mean([sess_avg{i}{1};sess_avg{i}{2};sess_avg{i}{3};sess_avg{i}{4};sess_avg{i}{5}])
         sub_avg_gsr{i}=mean([sess_avg_gsr{i}{1};sess_avg_gsr{i}{2};sess_avg_gsr{i}{3};sess_avg_gsr{i}{4};sess_avg_gsr`{i}{5}])

    end
end

allsub=[sub_avg{1};sub_avg{2};sub_avg{3};sub_avg{4};sub_avg{5};sub_avg{6};sub_avg{7};sub_avg{8};sub_avg{9};sub_avg{10};sub_avg{11}];
allsub_gsr=[sub_avg_gsr{1};sub_avg_gsr{2};sub_avg_gsr{3};sub_avg_gsr{4};sub_avg_gsr{5};sub_avg_gsr{6};sub_avg_gsr{7};sub_avg_gsr{8};sub_avg_gsr{9};sub_avg_gsr{10};sub_avg_gsr{11}];

avg_ICC=mean(allsub);
avg_ICC_GSR=mean(allsub_gsr);


avg_ICC_reshape=reshape(avg_ICC, 91,109,91);
avg_ICC_GSR_reshape=reshape(avg_ICC_GSR, 91,109,91);
diff_ICC_reshape=reshape(diff_ICC, 91,109,91);

save_avw(avg_ICC_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/AvgICC_withinsess_sub'), 'f', [2 2 2 2]);
save_avw(avg_ICC_GSR_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/AvgICC_GSR'), 'f', [2 2 2 2]);
save_avw(diff_ICC_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/Diff_AvgICC'), 'f', [2 2 2 2]);

%% compute average ICC across all subjects - concatenated scans.
all_ICC={};
all_ICC_GSR={};
%diff_ICC= [];
sess_avg=sessions;

for i = 1:11
    for j=1:nsess(i)
        n=numscans(j,i)
        all_ICC{j} = reshape(ICC_final_cat{i}{j}{n},[],1)';%one file per session 1,2,3,4 and 5 
        all_ICC_GSR{j} = reshape(ICC_GSR_final_cat{i}{j}{n},[],1)';
    end
    if i==6
         sub_avg{i}=mean([all_ICC{1};all_ICC{2};all_ICC{3};all_ICC{4}])
         sub_avg_GSR{i}=mean([all_ICC_GSR{1};all_ICC_GSR{2};all_ICC_GSR{3};all_ICC_GSR{4}])

    else 
         sub_avg{i}=mean([all_ICC{1};all_ICC{2};all_ICC{3};all_ICC{4};all_ICC{5}])
         sub_avg_GSR{i}=mean([all_ICC_GSR{1};all_ICC_GSR{2};all_ICC_GSR{3};all_ICC_GSR{4};all_ICC_GSR{5}])
    end
end

allsub_cat=[sub_avg{1};sub_avg{2};sub_avg{3};sub_avg{4};sub_avg{5};sub_avg{6};sub_avg{7};sub_avg{8};sub_avg{9};sub_avg{10};sub_avg{11}];
allsub_GSR_cat=[sub_avg_GSR{1};sub_avg_GSR{2};sub_avg_GSR{3};sub_avg_GSR{4};sub_avg_GSR{5};sub_avg_GSR{6};sub_avg_GSR{7};sub_avg_GSR{8};sub_avg_GSR{9};sub_avg_GSR{10};sub_avg_GSR{11}];

avg_ICC_cat=mean(allsub_cat);
avg_ICC_GSR_cat=mean(allsub_GSR_cat);

avg_ICC_cat_reshape=reshape(avg_ICC_cat, 91,109,91);
avg_ICC_GSR_cat_reshape=reshape(avg_ICC_GSR_cat, 91,109,91);
%diff_ICC_reshape=reshape(diff_ICC, 91,109,91);

save_avw(avg_ICC_cat_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/AvgICC_withinsess_sub_cat'), 'f', [2 2 2 2]);
save_avw(avg_ICC_GSR_cat_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/AvgICC_GSR_cat'), 'f', [2 2 2 2]);
%save_avw(diff_ICC_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/Diff_AvgICC_cat'), 'f', [2 2 2 2]);


%% Figures - averaged ICC between scans
figure(1) %heatmap without GSR
%handle_gm_GSR{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0,rand1000)');
caxis(handle_gm{i}{j}{k}, [-8 8])
handle_gm{i}{j}{k}.FontSize=6;
handle_gm{i}{j}{k}.Title='Timeseries values without GSR';
handle_gm{i}{j}{k}.GridVisible = 'off'
handle_gm{i}{j}{k}.Colormap = gray

figure(2) %heatmap with GSR
%handle_gm_GSR{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0,rand1000)');
caxis(handle_gm_GSR{i}{j}{k}, [-8 8])
handle_gm_GSR{i}{j}{k}.FontSize=6;
handle_gm_GSR{i}{j}{k}.Title='Timeseries values WITH GSR';
handle_gm_GSR{i}{j}{k}.GridVisible = 'off'
handle_gm_GSR{i}{j}{k}.Colormap = gray

figure(3) %plot correlation coefficients before and after GSR.
%Cupper_gm{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
ksdensity(Cupper_gsr{i}{j}{k})
hold on;
ksdensity(Cupper_gm{i}{j}{k})
title('GM Correlation coefficients before (red) and after (blue) GSR')

figure(4) %difference in voxel values before and after gm
%handle_gm_GSR_DIFF{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0, rand1000)'-ts_GM(sum(outlierframes,2)==0,rand1000)');
caxis(handle_gm_GSR_DIFF{i}{j}{k}, [-8 8])
handle_gm_GSR_DIFF{i}{j}{k}.FontSize=6;
handle_gm_GSR_DIFF{i}{j}{k}.Title='Difference in timeseries values with GSR';
handle_gm_GSR_DIFF{i}{j}{k}.GridVisible = 'off'
handle_gm_GSR_DIFF{i}{j}{k}.Colormap = gray

figure(5) %plot ICC values before and after GSR
ksdensity(ICC_GSR{i}{j}{k})
title('ICC values after GSR')
hold on;
ksdensity(ICC{i}{j}{k})
title('ICC values before (red) and after (blue) GSR')



%% Figures - concatenated
figure(6) %heatmap without GSR
%handle_gm_GSR{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0,rand1000)');
caxis(handle_gm_cat{i}{j}{k}, [-8 8])
handle_gm_cat{i}{j}{k}.FontSize=6;
handle_gm_cat{i}{j}{k}.Title='Timeseries values without GSR';
handle_gm_cat{i}{j}{k}.GridVisible = 'off'
handle_gm_cat{i}{j}{k}.Colormap = gray

figure(7) %heatmap with GSR
%handle_gm_GSR{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0,rand1000)');
caxis(handle_gm_GSR_cat{i}{j}{k}, [-8 8])
handle_gm_GSR_cat{i}{j}{k}.FontSize=6;
handle_gm_GSR_cat{i}{j}{k}.Title='Timeseries values WITH GSR';
handle_gm_GSR_cat{i}{j}{k}.GridVisible = 'off'
handle_gm_GSR_cat{i}{j}{k}.Colormap = gray

figure(8) %plot correlation coefficients before and after GSR.
%Cupper_gm{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
histogram(Cupper_gsr_cat{i}{j}{k})
hold on;
histogram(Cupper_gm_cat{i}{j}{k})
title('GM Correlation coefficients before (red) and after (blue) GSR - concatenated')

figure(9) %difference in voxel values before and after gm
%handle_gm_GSR_DIFF{i}{j}{k}=heatmap(ts_GSR(sum(outlierframes,2)==0, rand1000)'-ts_GM(sum(outlierframes,2)==0,rand1000)');
caxis(handle_gm_GSR_DIFF_cat{i}{j}{k}, [-8 8])
handle_gm_GSR_DIFF_cat{i}{j}{k}.FontSize=6;
handle_gm_GSR_DIFF_cat{i}{j}{k}.Title='Difference in timeseries values with GSR';
handle_gm_GSR_DIFF_cat{i}{j}{k}.GridVisible = 'off'
handle_gm_GSR_DIFF_cat{i}{j}{k}.Colormap = gray

figure(10) %plot ICC values before and after GSR
ksdensity(ICC_GSR_cat{i}{j}{k})
title('Concatenated - ICC values after GSR')
hold on;
ksdensity(ICC_cat{i}{j}{k})
title('ICC values before (red) and after (blue) GSR')



%% pairwise voxel correlation coefficients (GM) between concatenated scans and averaged scans.
figure(8) %plot correlation coefficients before and after GSR.
%Cupper_gm{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
i=1
j=1
%Cupper_gsr = with GSR
%Cupper_gm = without GSR
figure(1)
for i=1:4
    subplot(1,4,i)
    a=histogram(Cupper_gsr_cat{1}{1}{4}, 1000)
    a.FaceColor ='red'
    a.FaceAlpha = 0.2
    a.EdgeColor ='red'
    a.EdgeAlpha = 0.2
    a
    hold on;
    b=histogram(Cupper_gsr{1}{4}{i}, 1000)
    b.EdgeColor ='blue'
    b.EdgeAlpha = 0.2
    b.FaceColor ='blue'
    b.FaceAlpha = 0.2
    title(strcat('Concatenated vs Scan ', num2str(i)))
    b
end

set(gcf, 'Position', [10 10 2000 500])
sgtitle('Pairwise voxel timeseries correlation coefficients in subject 1 session 4')
legend('Concatenated', 'Not concatenated')


size(Cupbper_gsr{i}{j}{1})
size(Cupper_gsr_cat{i}{j}{4})

histogram(Cupper_gm_cat{i}{j}{4})
title('GM Correlation coefficients before (red) and after (blue) GSR - concatenated')

figure(3) %plot correlation coefficients before and after GSR.
%Cupper_gm{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
hold on;
ksdensity(Cupper_gm{i}{j}{k})



%% Boxplots of ICC of voxels connected vs unconnected to lesion area
traject = figure(1)
set(traject, 'Position', [0 0 1400 1200])
traject
mean_sess=[];
p =1;
for i=1:4
    if i ==6
        continue
    end
    
    for j=1:5
        scans=[];
        for k=1:numscans(j,i)
            scans=[scans; ICC{i}{j}{k}]; %masks out things not in the gray matter for plotting.
        end
        mean_sess{i,j}=mean(scans,1);
        pos=mean_sess{i,j}(lesion_dc{i}(logical(graymask{i}'),:)>3);
        sizepos=size(pos,2);
        zer=mean_sess{i,j}(lesion_dc{i}(logical(graymask{i}'),:)<3);
        sizezer=size(zer,2);
        grp=[zeros(1,size(pos,2)), ones(1,size(zer,2))];
        fun = [pos zer];
        hold on;
       subplot(5,4, p);
        boxplot(fun, grp);
        title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
        xlabel('Connectivity to lesion area')
        ylabel('ICC')
        xticklabels({'Connected', 'Not connected'})
        ylim([0 0.1])
        p=p+1;
    end
end     


traject = figure(2)
set(traject, 'Position', [0 0 1700 1200])
traject
mean_sess=[];
p =1;
for i =7:11
    for j=1:5
        scans=[];
        for k=1:numscans(j,i)
            scans=[scans; ICC{i}{j}{k}];
        end
        mean_sess{i,j}=mean(scans,1);
        pos=mean_sess{i,j}(lesion_dc{i}>3);
        sizepos=size(pos,2);
        zer=mean_sess{i,j}(lesion_dc{i}<3);
        sizezer=size(zer,2);
        grp=[zeros(1,size(pos,2)), ones(1,size(zer,2))];
        fun = [pos zer];
        hold on;
        subplot(2,5,p);
        boxplot(fun, grp);
        title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
        xlabel('Connectivity to lesion area')
        ylabel('ICC')
        xticklabels({'Connected', 'Not connected'})
        ylim([0 0.1])
        p=p+1;
    end
end    



%% correlation between ICC and strength of connectivity to lesion area, within areas structurally connected to the lesion.
% averaged across scans in sessions.
clear gm;
p=1;
for i=6
    subject = strcat('SUB', num2str(i));
    graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel

    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels>
    for j=1:nsess(i)
        scans=[];
        for k=1:numscans(j,i)
            scans=[scans; ICC_GSR{i}{j}{k}]; %masks out things not in the gray matter for plotting.
        end
       icc = mean(scans);
       lesion = lesion_dc{i}(graymask{i}==1);
       subplot(1,4,p);
       semilogy(icc(lesion>4), lesion(lesion>4),'b.');
       hold on;
       %plot(icc(lesion>1), log(lesion(lesion>1)), 'r.')
       %coeffs =polyfit(icc(lesion>5), log(lesion(lesion>5)), 1);
       %fitX = linspace(min(icc(lesion>5)), max(icc(lesion>5)), 200);
       %fitY = polyval(coeffs, fitX);
       %hold on; 
    %plot(fitX, fitY, 'b-', 'LineWidth', 3);
    title(['subject ', num2str(i), 'session ', num2str(j)])
    sgtitle('with GSR. threshold = 4. scans are averaged within sessions')
    xlabel('ICC')
    ylabel('connectivity to lesion (log)')
    xlim([0.02 0.07])
    p=p+1;
    end
end
set(gcf, 'Position', [10 10 1200 1200])

export_fig(strcat('/mnt/shared_data3/emo4002/pons_sfmodelling/figures/231119_correlation_ICC_disconnectivity/sub1_5_withinarea_correlation.pdf'))


%% relationship between lesionTract threshold & correlation with ICC.

%calculate mean correlation w/ ICC at different thresholds for connectivity
%with ICC (lesion>X)
clear gm;
p=1;
corr_i=[]
corr_subj =[];
for i=1:11
    subject = strcat('SUB', num2str(i));
    graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel

    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels
    corr_sess=[];
    for j=1:nsess(i) % num sessions per subject = always 5 except subject 6 has 0.
        scans=[];
        corr_k=[];
        for k=1:numscans(j,i) % num scans per session = usually 3, sometimes 4 
            corrx=[];
            scans=[scans; ICC{i}{j}{k}]; %masks out things not in the gray matter for plotting.
            icc = scans(k,:)';
            lesion = lesion_dc{i}(graymask{i}==1);
            for x = 1:50
                corrx(x) = corr(icc(lesion>x), log(lesion(lesion>x))); %corrx_sp = spearman

            end
            corr_k = vertcat(corr_k, corrx); 
        end
        corr_j = mean(corr_k,1); %mean correlation across scans within session j.
        corr_sess = vertcat(corr_sess, corr_j);
    end
    corr_i=mean(corr_sess,1)%mean across sessions within subject i.
    corr_subj =vertcat(corr_subj, corr_i)
end

corrICC_disc = mean(corr_subj,1)

%plot
mean_corr=corr_subj
plot(mean_corr', 'r.')
title({'correlation between ICC and structural connectivity to ','lesion area with different lesion connectivity thresholds'})
xlabel('threshold for connectivity')
ylabel('correlation (pearsons)')



%% plot number of voxels included in lesion area at different thresholds for connectivity
%with ICC (lesion>X)
clear gm;
p=1;
corr_i=[]
corr_subj =[];
for i=1:11
    subject = strcat('SUB', num2str(i));
    graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/avg/wc1c',subject,'_Savg_GM_binarized.nii.gz'));
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel

    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels
    corr_sess=[];
    for j=1 % num sessions per subject = always 5 except subject 6 has 0.
        scans=[];
        corr_k=[];
        for k=1 % num scans per session = usually 3, sometimes 4 
            corrx=[];
            scans=[scans; ICC{i}{j}{k}]; %masks out things not in the gray matter for plotting.
            icc = scans(k,:)';
            lesion = lesion_dc{i}(graymask{i}==1);
            for x = 1:50
                corrx(x) = size(lesion(lesion>x), 1); %corrx_sp = spearman
            end
            corr_k = vertcat(corr_k, corrx); 
        end
        corr_sess = vertcat(corr_sess, corr_k);
    end
    corr_subj =vertcat(corr_subj, corr_sess)
end

%plot
mean_corr=corr_subj
plot(mean_corr', 'r.')
title({'number of voxels "connected" to ','lesion area with different lesion connectivity thresholds'})
xlabel('threshold for connectivity')
ylabel('number of voxels')


%% generate gray matter control mask.
for i=19:44
    if i==36
        continue
    end
    if i==39
        continue
    end
    subject = strcat('SUB', num2str(i));
    graymask{i} = read_avw(strcat('/home/emo4002/colossus_shared3/control_subs/',subject,'/avg/wc1c',subject,'_Savg.nii'));
    graymask{i}=reshape(graymask{i},[],1);
    graymask{i} = graymask{i}>0;
end

intersection=graymask{19}.*graymask{20}.*graymask{21}.*graymask{22}.*graymask{23}.*graymask{24}.*graymask{25}.*graymask{26}.*graymask{27}.*graymask{28}.*graymask{29}.*graymask{30}.* ...
    graymask{31}.*graymask{32}.*graymask{33}.*graymask{34}.*graymask{35}.*graymask{37}.*graymask{38}.*graymask{40}.*graymask{41}.*graymask{42}.*graymask{43}.*graymask{44};

mata=[graymask{19},graymask{20},graymask{21},graymask{22},graymask{23},graymask{24},graymask{25},graymask{26},graymask{27},graymask{28},graymask{29},graymask{30}, ...
  graymask{31},graymask{32},graymask{33},graymask{34},graymask{35},graymask{37},graymask{38},graymask{40},graymask{41},graymask{42},graymask{43},graymask{44}];

p=0;
for i=1:size(mata,1)
    if sum(mata(i,:),2)==23
        mata(i,:)=boolean([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
        p=p+1;
    end
    if sum(mata(i,:),2)==22
        mata(i,:)=boolean([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
        p=p+1;
    end
    if sum(mata(i,:),2)==21
        mata(i,:)=boolean([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
        p=p+1;
    end
    if sum(mata(i,:),2)==20
        mata(i,:)=boolean([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]);
        p=p+1;
    end
end

imagesc(mata)
hist(sum(mata,2),1000)
hold on;
ylim([0 40000]);
xlabel('Number of subjects included');
ylabel('Number of subjects with data at a voxel')

int_23=prod(mata,2);
int_22=prod(mata,2);

intersect = reshape(intersection,91,109,91);
intersect_22 = reshape(int_22,91,109,91);
save_avw(intersect,strcat('/home/emo4002/colossus_shared3/control_subs/control_intersection'),'f',[2 2 2 3]);
save_avw(intersect_22,strcat('/home/emo4002/colossus_shared3/control_subs/control_intersection_20_24'),'f',[2 2 2 3]);




%% dump

Vgray=read_avw('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/avg/wc1cSUB1_Savg_GM_binarized.nii.gz'); %3D volume
Vgray1=reshape(Vgray,[],1); %flattened 1D vector with a 1/0 for each voxel

Vlesion=read_avw('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/lesionTract/HCP_MGH_32fold/wc_SUB1_lesionTract_GM_masked.nii.gz');
Vlesion=reshape(Vlesion,[],size(Vlesion,4)); %flattened 2D matrix that is <voxels>
Vlesion=Vlesion(logical(Vgray),:); %same but now <grayvoxels>

Vts=read_avw('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/func/S1/denoise_swaufunc1_GM_masked.nii.gz');
Vts1=reshape(Vts,[],size(Vts,4)); %flattened 2D matrix that is <voxels>

Vts1_gray=Vts1(logical(Vgray1),:); %same but now <grayvoxels> x <time>
ts=Vts1_gray.'; %flip to <time> x <grayvoxels>

outliers=load('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/func/S1/art_regression_outliers_aufunc_1.mat');
outlierframes = outliers.mat;
outlierframes = outlierframes(:, sum(outlierframes,1) > 0); %
%find(sum(outlierframes, 2))

motionframes =find(sum(outlierframes,2));
motionframes=unique([motionframes; motionframes+1]);

outlierframes=zeros(124, size(motionframes,1));

%create new outlierframes matrix adding frame before and after frame
%identified by CONN
for i=1:size(motionframes,1)
    outlierframes(motionframes(i),i) = 1;
end

% confound regression
meants=mean(ts,2); %mean across all gray matter voxels
confounds_GSR=[outlierframes meants [0; diff(meants)]];
Q_GSR=eye(size(confounds_GSR,1))-confounds_GSR*pinv(confounds_GSR);
ts_GSR=Q_GSR*ts;

ts_clean_GSR=ts_GSR(~sum(outlierframes,2),:); %only timeseries without motion
ts_clean=ts(~sum(outlierframes,2),:);

ICC_GSR =intrinsic_connectivity_contrast(ts_clean_GSR);
ICC = intrinsic_connectivity_contrast(ts_clean);

Vlesion=Vlesion';
plot(ICC, Vlesion, 'o');



%add FSL libraries
if(isempty(which('save_avw')))
        addpath([getenv('FSLDIR') '/etc/matlab']);
    end
libPath =getenv('LD_LIBRARY_PATH'); libPath =[libPath pathsep '/usr/lib/fsl/5.0']; setenv('LD_LIBRARY_PATH',libPath);

%load GM timeseries.
epi=read_avw('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/func/S1/denoise_swaufunc1_GM_masked.nii.gz');
%load outlier mat file generated by conn
outliers=load('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/func/S1/art_regression_outliers_aufunc_1.mat')

%remove outliers from GM timeseries
outlier_vector=sum(outliers.mat,2); %new vector; outlier rows =1, everything else =0
outlier_log = logical(outlier_vector); %turn into logical to do logical indexing
epi=epi(:, :, :, ~outlier_log); %only keep volumes that are not outliers.
size(epi)%confirm reduced # volumes.

save_avw(epi, '/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/func/S1/denoise_swaufunc1_GM_masked_outliers.nii.gz','f',[2 2 2 2])
clear all;
conn_vol('/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/func/S1/denoise_swaufunc1_GM_masked.nii.gz'); %this is a Nt x (maskvoxels)

save('/home/emo4002/colossus_shared3/pons_sfmodelling/epi.mat')

test =intrinsic_connectivity_contrast(epi)

fsl_glm
path_tract='/home/emo4002/colossus_shared3/pons_sfmodelling/SUB1/lesionTract/HCP_MGH_32fold/wc_GM_SUB1_lesionTract.nii.gz';
%should really be S1 not SUB1.
matrix_tract = voxel_intensity_from_nii(path_tract);

path_icc='/home/emo4002/colossus_shared3/pons_sfmodelling/conn_projects/conn_project_single_subejct1_denoised/results/firstlevel/V2V_01/wc_ICC_GM_SUB1.nii.gz';
matrix_icc = voxel_intensity_from_nii(path_icc);

corrcoef(matrix_tract, matrix_icc)
Vdata3,outfile,'f',voxdim
