function [] = calculateZScores(nscans, nsess, studydir, resultsdir);
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

    load(strcat(studydir, resultsdir, 'ICC_GSR_23subjects.mat'),'ICC_GSR_cat')
    load(strcat(studydir, resultsdir, 'crtl_cat_GSR_mean.mat'), 'ctrl_cat_GSR_mean')
    load(strcat(studydir, resultsdir, 'crtl_cat_GSR_stdev.mat'), 'ctrl_cat_GSR_stdev')

    SUBzscore=[];
    for i=1:23
        for j=1:nsess(i)        
            % load mean and stddev of controls
            SUBi_ctrlmean = ctrl_cat_GSR_mean; 
            SUBi_ctrlstdev = ctrl_cat_GSR_stdev;
        
            % subtract control mean from subject's ICC, divide by sd of controls
            %z-score calculation: (Xstrokept-mean of controls)/(standard deviation of ctls)
            SUBzscore{i,j}=(ICC_GSR_cat{i}-SUBi_ctrlmean)./SUBi_ctrlstdev
            GM = read_avw(str_cat(studydir, 'c1referenceT1.nii'));
            GM_reshape = reshape(GM, [1 902629])

            % reshape data to 3D
            p=1; %counter
            zICC_GSR=[];
            zICC_GSR_3D=[];
            for z=1:size(GM_reshape,1) %starts with the gray matter mask(which has all voxels in the 3d volume.)
                % makes gray matter voxels equal to their ICC value.
                if GM_reshape(z)==1
                    zICC_GSR(z)=ICC_GSR_cat{i}{j}(p); 
                    p=p+1;
                else
                    zICC_GSR(z) =0;
                end
            end
            zICC_GSR_3D = reshape(zICC_GSR,91,109,91);
            save_avw(zICC_GSR_3D,strcat(studydir, resultsdir, 'SUB',num2str(i), '_S',num2str(j), '_zscoreICC'),'f',[2 2 2 2]);

        end
    end
end


%% Boxplots of ICC voxels connected vs unconnected to lesion area - concatenated scans, Z-SCORED!
traject = figure(1)
set(traject, 'Position', [0 0 1400 1200])
traject
p =1;
for i=1:3
    subject = strcat('SUB',num2str(i));
    %load gray matter masks
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5 %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>

    %lesion_dc{i}=lesion_dc{i}.*logical(graymask{i}); %same but now <grayvoxels>
    for j=1:nsess(i)
        pos=SUBzscore{i,j}(lesion_dc{i}(graymask{i})>1);
        sizepos=size(pos,2);
        zer=SUBzscore{i,j}(lesion_dc{i}(graymask{i})<1);
        sizezer=size(zer,2);
        grp=[zeros(1,size(pos,2)), ones(1,size(zer,2))];
        fun = [pos zer];
        hold on;
        subplot(4,5, p);
        boxplot(fun, grp);
        ylim([-6 6])
        title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
        xlabel('Connectivity to lesion area')
        ylabel('z-score')
        xticklabels({'Connected', 'Not connected'})
        p=p+1;
    end
end


histogram(SUBzscore{i,j})
xlim([-5 5])

figure(2)
connected=SUBzscore{i,j}(lesion_dc{i}(graymask{i}==1)>0);
disconnected = SUBzscore{i,j}(lesion_dc{i}(graymask{i}==1)==0);
histogram(connected,1000, 'FaceColor', 'r', 'EdgeAlpha', 0.2)
hold on;
histogram(disconnected,100000, 'FaceColor', 'b', 'EdgeAlpha', 0.2)
xlim([-5 5])
violinplot(connected, disconnected)

figure(3)
histogram(SUBzscore{i,j}(lesion_dc{i}(graymask{i}==1)>5),1000, 'FaceColor', 'r', 'EdgeAlpha', 0.2)
xlim([-5 5])
hold on;
histogram(SUBzscore{i,j}(lesion_dc{i}(graymask{i}==1)<5),10000, 'FaceColor', 'g', 'EdgeAlpha', 0.2)
xlim([-5 5])
title('Histogram of ICC z-scores in brain areas connected (red) vs unconnected (green) to lesion. ICCs of zero not yet fixed');

%% boxplots of ICC connected to lesion mask session 1 vs session 5 and not connected session 1 to ssession 5

traject = figure(6)
set(traject, 'Position', [0 0 2000 1200])
traject
d=1
for i=1:11
    subject = strcat('SUB',num2str(i));
    %load gray matter masks
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5 %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    baseline_z = SUBzscore{i,1}(lesion_dc{i}(graymask{i})<1);% connected
    if i==6
        final_z = SUBzscore{i,4}(lesion_dc{i}(graymask{i})<1);% connected
    else
        final_z = SUBzscore{i,5}(lesion_dc{i}(graymask{i})<1);% connected
    end
    [h,p,ci,stats] = ttest2(baseline_z, final_z);
    tstat=stats.tstat
    tst(i)=tstat
    p_n(i)=p
    [h,p,ci,stats]=ttest2(baseline_z,final_z,'Vartype','unequal') %asssume unequal variance.
    tst_w(i)=stats.tstat
    p_w(i)=p
    subplot(2,11,d)
    pl=[baseline_z',final_z'];
    boxplot(pl)
    d=d+1;
    title(strcat('t-stat=',num2str(tst(i))))
    hold on;
    %significance stuff
    % yt = get(gca, 'YTick');
    %   axis([xlim    0  ceil(max(yt)*1.2)])
    %  xt = get(gca, 'XTick');
    hold on
    %   plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
    %    hold off
    ylim([-20 20])
end
%
for i=1:11
    subject = strcat('SUB',num2str(i));
    %load gray matter masks
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5 %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    baseline_z = SUBzscore{i,1}(lesion_dc{i}(graymask{i})>1);% connected
    if i==6
         final_z = SUBzscore{i,4}(lesion_dc{i}(graymask{i})>1);% connected
    else
         final_z = SUBzscore{i,5}(lesion_dc{i}(graymask{i})>1);% connected
    end
    [h,p,ci,stats] = ttest2(baseline_z, final_z);
    tstat=stats.tstat
    tst_conn(i)=tstat
    p_n(i)=p
    [h,p,ci,stats]=ttest2(baseline_z,final_z,'Vartype','unequal') %asssume unequal variance.
    tst_w(i)=stats.tstat
    p_w(i)=p    
    subplot(2,11,d)
    pl=[baseline_z',final_z'];
    boxplot(pl)
    d=d+1;
    title(strcat('t-stat=',num2str(tst(i))))
    hold on;
    %significance stuff
   % yt = get(gca, 'YTick');
 %   axis([xlim    0  ceil(max(yt)*1.2)])
  %  xt = get(gca, 'XTick');
    hold on
 %   plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
%    hold off
ylim([-20 20])
end

d=1
%just connected regionss
for i=1:11
    subject = strcat('SUB',num2str(i));
    %load gray matter masks
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5 %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    baseline_z = SUBzscore{i,1}(lesion_dc{i}(graymask{i})>1);% connected
    if i==6
        final_z = SUBzscore{i,4}(lesion_dc{i}(graymask{i})>1);% connected
    else
        final_z = SUBzscore{i,5}(lesion_dc{i}(graymask{i})>1);% connected
    end
    [h,p,ci,stats] = ttest2(baseline_z, final_z);
    tstat=stats.tstat
    tst_conn(i)=tstat
    p_n(i)=p
    [h,p,ci,stats]=ttest2(baseline_z,final_z,'Vartype','unequal') %asssume unequal variance.
    tst_w(i)=stats.tstat
    p_w(i)=p
    subplot(1,11,d)
    pl=[baseline_z',final_z'];
    boxplot(pl)
    d=d+1;
    title(strcat('t-stat=',num2str(tst_conn(i))))
    hold on;
    %significance stuff
    % yt = get(gca, 'YTick');
    %   axis([xlim    0  ceil(max(yt)*1.2)])
    %  xt = get(gca, 'XTick');
    hold on
    %   plot(xt([1 2]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 2])), max(yt)*1.15, '*k')
    %    hold off
  %  ylim([-5 10])
    xticklabels({'S1','S5'})
   % ylabel('z-scored ICC in areas connected to lesion')
end


scatter(tst,tst_conn)
%% t-test between ICC of voxels connected to the lesion mask on session 1 vs session 5
% & correlation with recovery scores

clf;
tst=[];
rec=[];

traject = figure(4)
set(traject, 'Position', [0 0 1200 1200])
traject;

for i=1:11
    subject = strcat('SUB',num2str(i));
    %load gray matter masks
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5; %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    df =readmatrix('/mnt/shared_data3/emo4002/pons_sfmodelling/stroke_pts/demographics/Demographics_11pts_fuglmeyer_avg.csv');
    baseline_z = SUBzscore{i,1}(lesion_dc{i}(graymask{i})>1);% connected
    if i==6
         final_z = SUBzscore{i,4}(lesion_dc{i}(graymask{i})>1);% connected
    else
         final_z = SUBzscore{i,5}(lesion_dc{i}(graymask{i})>1);% connected
    end
    
  %  histogram(baseline_z)
  %  hold on;
  %  histogram(final_z)
    [h,p,ci,stats] = ttest2(final_z,baseline_z);
    tstat=stats.tstat
    tst(i)=tstat
    p_n(i)=p
    [h,p,ci,stats]=ttest2(baseline_z,final_z,'Vartype','unequal')
    tst_w(i)=stats.tstat
    p_w(i)=p
    if i==6
      recovery = df(i+1,5)-df(i+1,2)
    else
      recovery = df(i+1,6)-df(i+1,2)
    end
    rec(i)=recovery
end

plot(tst,rec, '.r', 'MarkerSize', 18)
hold on;
b=polyfit(tst,rec,1);
a=polyval(b,tst)
plot(tst,a)

[rho,p]=corr(tst', rec', 'Type', 'Pearson')
pears_rho = rho
pears_p=p
[rho,p]=corr(tst', rec', 'Type', 'Spearman')
spear_rho=rho
spear_p=p

xlabel('T-statistic: ICC Session 5 vs. ICC Sesssion 1', 'fontsize', 18)
ylabel('Session 5 - Session 1 Fugl-Meyer score', 'fontsize', 18)

ax =gca
ax.FontSize=18
txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
txt2=({strcat('Spearmans Rank Correlation = ', num2str(spear_rho)), strcat('p = ',num2str(round(spear_p,4)))})
text(-20,70,txt, 'fontsize', 18)
%text(-20,60,txt2, 'fontsize', 18)
title({'Relationship between motor recovery and change in session 1 vs session 5 ICC',' in cortical areas structurally connected to lesion'}, 'fontsize', 18)
saveas(gcf,'/home/emo4002/colossus_shared3/pons_sfmodelling/figures/181219_motorrecovery_vs_ICCchange/motorrecov_changeinICC_sess1vsess5_11subjects.png')

%% boxplot of S5-S1, S4-S1, S3-S1, etc. and motor recovery
rec=[];
tst=[];
clf;

traject = figure(5)
set(traject, 'Position', [0 0 1000 1000])
traject;

for i=1:11
    figure(i)
    subject = strcat('SUB',num2str(i));
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5; %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    df =readmatrix('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_fuglmeyer_avg.csv');
    if i==6
         z_score_dis=[SUBzscore{i,1}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,2}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,3}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,4}(lesion_dc{i}(graymask{i})>1)']
         boxplot(z_score_dis)
         xticklabels({'S1','S2','S3','S4'})
         title(strcat('Subject ', num2str(i)))
         xlabel('Sesssion')
         ylabel('Z-scored ICC values')
    else 
         z_score_dis=[SUBzscore{i,1}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,2}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,3}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,4}(lesion_dc{i}(graymask{i})>1)', SUBzscore{i,5}(lesion_dc{i}(graymask{i})>1)']
         boxplot(z_score_dis)
         xticklabels({'S1','S2','S3','S4', 'S5'})
          title(strcat('Subject ', num2str(i)))
         xlabel('Sesssion')
         ylabel('Z-scored ICC values')
    end
    saveas(gcf,strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/figures/171219_icc_trajectory/sub',num2str(i),'.png'))
end


%% t-test between ICC of voxels connected to the lesion mask & correlation with recovery scores between SUBSESQUENT timepoints
clf;
tst=[];
rec=[];

traject = figure(5)
set(traject, 'Position', [0 0 1200 1200])
traject;
t=1;
for w=1:11
    for r=1:(nsess(w)-1)
        subject = strcat('SUB',num2str(w));
        graymask{w} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
        graymask{w}=reshape(graymask{w},[],1); %flattened 1D vector with a 1/0 for each voxel
        graymask{w}=graymask{w}>0.5; %load lesion tracts
        lesion_dc{w} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
        lesion_dc{w}=reshape(lesion_dc{w},[],size(lesion_dc{w},4)); %flattened 2D matrix that is <voxels>
        df =readmatrix('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_fuglmeyer_avg.csv');
        Sz_base= SUBzscore{w,r}(lesion_dc{w}(graymask{w})>1);% connected
        Sz_follow=SUBzscore{w,r+1}(lesion_dc{w}(graymask{w})>1);% connected
        [h,p,ci,stats] = ttest2(Sz_follow,Sz_base);
        tst(t)=stats.tstat
        recovery = df(w+1,r+2)-df(w+1,r+1)
        rec(t)=recovery
        t=t+1;
    end
   
  %  figure(w)
  %  plot(tst,rec, 'or')
  %  hold on;
  %  b=polyfit(tst,rec,1);
  %  a=polyval(b,tst)
  %  plot(tst,a);
  %  xlabel('t-statistic of ICC value distributions between sessions')
  %  ylabel('Between-session change in motor score')
  %  title(strcat('Subject ',num2str(w)))
  % saveas(gcf,strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/figures/171219_icc_trajectory/corr_motor_deltaICC_sub',num2str(w),'.png'));
end

plot(tst,rec, '.r', 'MarkerSize', 18) 
hold on;
b=polyfit(tst,rec,1);
a=polyval(b,tst)
plot(tst,a);

txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
%txt2=({strcat('Spearmans Rank Correlation = ', num2str(spear_rho)), strcat('p = ',num2str(round(spear_p,4)))})
text(-20,70,txt, 'fontsize', 18)
%text(-20,60,txt2, 'fontsize', 18)
[rho,p]=corr(tst', rec', 'Type', 'Pearson')
pears_rho = rho
pears_p=p
ax =gca
ax.FontSize=18
%[rho,p]=corr(tst', rec', 'Type', 'Spearman')
%%spear_rho=rho
%spear_p=p
xlabel('T-statistic: ICC between sequential sessions (Follow-up vs. Baseline)', 'fontsize', 18)
ylabel('Sequential session change in Fugl-Meyer score (Follow-up - Baseline)', 'fontsize', 18)
title({'Relationship between motor recovery & change in sequential baseline v. followup','ICC in cortical areas structurally connected to lesion'}, 'fontsize', 18)
saveas(gcf,'/home/emo4002/colossus_shared3/pons_sfmodelling/figures/181219_motorrecovery_vs_ICCchange/motorrecov_subsequentsess_changeinICC_11subjects.png')

%%

clf;
figure(4)
plot(tst,rec, 'or')
hold on;
b=polyfit(tst,rec,1);
a=polyval(b,tst)
plot(tst,a)

for r=1:4
    i=10; %SUBJECT 10
    subject = strcat('SUB',num2str(i));
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5; %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    df =readmatrix('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_fuglmeyer_avg.csv');
    Sz_base= SUBzscore{i,r}(lesion_dc{i}(graymask{i})>1);% connected
    Sz_follow=SUBzscore{i,r+1}(lesion_dc{i}(graymask{i})>1);% connected
    [h,p,ci,stats] = ttest2(Sz_base, Sz_follow);
    tst(r)=stats.tstat
    recovery = df(i+1,r+2)-df(i+1,r+1)
    Sz{r}=Sz_base
    rec(r)=recovery
end

clf;
figure(4)
plot(tst,rec, 'or')
hold on;
b=polyfit(tst,rec,1);
a=polyval(b,tst)
plot(tst,a)

%% correlation betwen change of ICC basesline vs followup (z-scores) vs structural connectivity to lesion area.
d=1
for i=1:5
     subject = strcat('SUB',num2str(i));
    %load gray matter masks
    graymask{i} = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii');
    graymask{i}=reshape(graymask{i},[],1); %flattened 1D vector with a 1/0 for each voxel
    graymask{i}=graymask{i}>0.5; %load lesion tracts
    lesion_dc{i} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
    lesion_dc{i}=reshape(lesion_dc{i},[],size(lesion_dc{i},4)); %flattened 2D matrix that is <voxels>
    lesion_dc{i}=lesion_dc{i}(graymask{i}>0.5)
    baseline_z = SUBzscore{i,1};% connected
    subplot(1,4,d)
    if i==6
         final_z = SUBzscore{i,4};% connected
    else
         final_z = SUBzscore{i,5};% connected
    end
    diff_ICC =baseline_z-final_z;
    subplot(1,4,d)
    histogram(diff_ICC(lesion_dc{i}>1), 'EdgeAlpha', 0.2, 'EdgeColor', 'red')
    hold on;
    plot([0, 0], [0,100]);
    hold on;
    yyaxis right
    histogram(diff_ICC(lesion_dc{i}<1), 'EdgeAlpha', 0.2, 'EdgeColor', 'blue')
    d=d+1;
end


%% Boxplots of ICC of voxels connected vs unconnected to lesion area -averaged scans, not z-scored
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
     %   all_ICC{j} = reshape(ICC_final_cat{i}{j}{n},[],1)';%one file per session 1,2,3,4 and 5 
        all_ICC_GSR{j} = reshape(ICC_GSR_final_cat{i}{j}{n},[],1)';
    end
    if i==6
    %     sub_avg{i}=mean([all_ICC{1};all_ICC{2};all_ICC{3};all_ICC{4}])
         sub_avg_GSR{i}=mean([all_ICC_GSR{1};all_ICC_GSR{2};all_ICC_GSR{3};all_ICC_GSR{4}])

    else 
         %sub_avg{i}=mean([all_ICC{1};all_ICC{2};all_ICC{3};all_ICC{4};all_ICC{5}])
         sub_avg_GSR{i}=mean([all_ICC_GSR{1};all_ICC_GSR{2};all_ICC_GSR{3};all_ICC_GSR{4};all_ICC_GSR{5}])
    end
end

%allsub_cat=[sub_avg{1};sub_avg{2};sub_avg{3};sub_avg{4};sub_avg{5};sub_avg{6};sub_avg{7};sub_avg{8};sub_avg{9};sub_avg{10};sub_avg{11}];
allsub_GSR_cat=[sub_avg_GSR{1};sub_avg_GSR{2};sub_avg_GSR{3};sub_avg_GSR{4};sub_avg_GSR{5};sub_avg_GSR{6};sub_avg_GSR{7};sub_avg_GSR{8};sub_avg_GSR{9};sub_avg_GSR{10};sub_avg_GSR{11}];

%avg_ICC_cat=mean(allsub_cat);
avg_ICC_GSR_cat=mean(allsub_GSR_cat);

%avg_ICC_cat_reshape=reshape(avg_ICC_cat, 91,109,91);
avg_ICC_GSR_cat_reshape=reshape(avg_ICC_GSR_cat, 91,109,91);
%diff_ICC_reshape=reshape(diff_ICC, 91,109,91);

%save_avw(avg_ICC_cat_reshape, strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/AvgICC_withinsess_sub_cat'), 'f', [2 2 2 2]);
save_avw(avg_ICC_GSR_cat_reshape, strcat('/home/emo4002/colossus_shared3/stroke_pts/AvgICC_GSR_cat'), 'f', [2 2 2 3]);
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


%Cupper_gsr = with GSR
%Cupper_gm = without GSR
figure(1)
for i=1:3
    subplot(1,3,i)
    a=histogram(Cupper_gsr_cat{2}{1}{3}, 1000)
    a.FaceColor ='red'
    a.FaceAlpha = 0.2
    a.EdgeColor ='red'
    a.EdgeAlpha = 0.2
    a
    hold on;
    b=histogram(Cupper_gsr{2}{3}{i}, 1000)
    b.EdgeColor ='blue'
    b.EdgeAlpha = 0.2
    b.FaceColor ='blue'
    b.FaceAlpha = 0.2
    title(strcat('Concatenated vs Scan ', num2str(i)))
    b
end

set(gcf, 'Position', [10 10 2000 500])
sgtitle('Pairwise voxel timeseries correlation coefficients in subject 2 session 1')
legend('Concatenated', 'Not concatenated')


size(Cupper_gsr{i}{j}{1})
size(Cupper_gsr_cat{i}{j}{4})

histogram(Cupper_gm_cat{i}{j}{4})
title('GM Correlation coefficients before (red) and after (blue) GSR - concatenated')

figure(3) %plot correlation coefficients before and after GSR.
%Cupper_gm{i}{j}{k}=C_gm(trimask); %should be %1000*1000*/2-1000
hold on;
ksdensity(Cupper_gm{i}{j}{k})


%% pairwise voxel correlation coefficients between concatenating and averaging 
figure(9)
a=histogram(Cupper_gsr_cat{1}{3}{3}, 1000)
a.FaceColor ='red'
a.FaceAlpha = 0.2
a.EdgeColor ='red'
a.EdgeAlpha = 0.2
a
hold on;
b=histogram(mean([Cupper_gsr{1}{3}{1},Cupper_gsr{1}{3}{2},Cupper_gsr{1}{3}{3}], 2), 1000)
b.EdgeColor ='blue'
b.EdgeAlpha = 0.2
b.FaceColor ='blue'
b.FaceAlpha = 0.2
b;
title(strcat('Concatenated vs Scan ', num2str(i)))

set(gcf, 'Position', [10 10 500 500])
legend('Concatenated', 'Averaged')




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
