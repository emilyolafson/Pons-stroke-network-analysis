function [tst, rec, nvoxels] = makeFigures(nsess, studydir, figuresdir, resultsdir, disconnectivitydir, figs)
% makeFigures: generate figures for Pontine stroke network analysis project. 
% INPUT: 
%     nsess: number of sessions per subject (columns)
%     studydir: folder containing subject folders and folder for outputs.
%     figures: name of folder where figures will be saved.
%     results: name of folder where ICC values are.
%     disconnectivitydir: name of folder where disconnectivity maps are.
%     figs: list of figures to be generated.
%
    % load z-scores
    load(strcat(studydir, resultsdir, 'zICC_23subjects.mat'),'SUBzscore')
    
    % load gray matter masks
    GM = read_avw(strcat(studydir, 'c1referenceT1.nii')); 
    GM_reshape = reshape(GM, [902629 1]);
    GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
    GM_reshape(GM_reshape <= 0.25) = 0;
    GM_reshape = logical(GM_reshape);
    df = readmatrix(strcat(studydir, 'demog_strokepts2.csv'));
    masks=[{'cortex'}, {'brainstem'}, {'cerebellum'}];

    % set analysis variables
    thresh = 0.5; % value above which voxels are considered to be 'connected' to the lesion

    % Figure 1 - boxplots of voxels connected vs unconnected to lesion
    if ismember(1, figs)
        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1400 280])

        for i = 1:23
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i}{j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i}{j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(1,nsess(i), j);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
            end
            saveas(gcf, strcat(studydir, figuresdir, 'boxplots_connect-vs-disconnect_', num2str(i), '.png'));
        end
    end
    close all;


    % Figure 2 - change in motor scores vs change in ICC from the last scan to the first scan 
    if  ismember(2, figs) 
        tst=[]; %tstatistic (change in ICC)
        rec=[]; %recovery (change in fugl-meyer)

        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1000 1000]);

        for i=1:23
            subject = strcat('SUB',num2str(i));

            %load disconnectivity files
            lesion_dc=read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
                
            baseline_z = SUBzscore{i}{1}(lesion_dc(GM_reshape)>thresh);

            if i==6
                final_z = SUBzscore{i}{4}(lesion_dc(GM_reshape)>thresh);
            elseif i==12
                final_z = SUBzscore{i}{3}(lesion_dc(GM_reshape)>thresh);
            elseif i==20
                final_z = SUBzscore{i}{2}(lesion_dc(GM_reshape)>thresh);
            else % most subjects came for 5 sessions.
                final_z = SUBzscore{i}{5}(lesion_dc(GM_reshape)>thresh);
            end
            
            % calculate the difference in z-score normalized ICC values at baseline and followup.
            % metric is a t-statistic using a paired t-test with unequal variances.
            % positive t-stat means distribution 
            [h,p,ci,stats]=ttest2(baseline_z,final_z,'Vartype','unequal')
            tst(i)=stats.tstat
            p(i)=p

            % load motor recovery scores (final - initial)
            % column 4 = sesssion 1
            % column 8 = sesssion 5
            if i==6 % only 4 sessions (sessions 1, 2, 3, and 4)
                recovery = df(i,7)-df(i,4)
            elseif i==12 % only 3 sessions (session 1, 2, and 3)
                recovery = df(i,6)-df(i,4)
            elseif i==20 % only 2 sessions (session 1 and 2)
                recovery = df(i,5)-df(i,4)
            else
                recovery = df(i,8)-df(i,4)
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
        ax = gca
        ax.FontSize = 18
        txt = ({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        %txt2 = ({strcat('Spearmans Rank Correlation = ', num2str(spear_rho)), strcat('p = ',num2str(round(spear_p,4)))})
        text(-20,70,txt, 'fontsize', 18)
        title({'Relationship between motor recovery and change in session 1 vs session 5 ICC',' in cortical areas structurally connected to lesion'}, 'fontsize', 18)
        saveas(gcf, strcat(studydir, figuresdir, '/correlation_allGM_changeICCvschangeFuglMeyer_baseline-vs-lastFU.png'))
    end
    close all;

    % Figure 3 - Correlation betwen inter-session change in ICC &inter-session change in Fugl-Meyer
    % also can be region-specific
    if ismember(3, figs)
        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1000 1000])

        tst=[]; %tstatistic (change in ICC)
        rec=[]; %recovery (change in fugl-meyer)

        t=1;
        for i=1:23
            if i==13
                continue
            end
            disp(num2str(i))
            %load disconnectivity files
            lesion_dc=read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            k=3;
            mask = getMask(masks{k});
            mask = reshape(mask, [902629 1]);
            intersect = GM_reshape.*logical(mask);
            lesion_intersect = lesion_dc.*intersect;
              
            if (i==6)
                for r=1:3
                    baseline_z = SUBzscore{i}{r}(lesion_intersect(GM_reshape==1)>thresh);
                    final_z = SUBzscore{i}{r+1}(lesion_intersect(GM_reshape==1)>thresh);
                    recovery = df(i,r+4)-df(i,r+3);
                    [h,p,ci,stats] = ttest2(baseline_z,final_z);
                    tst(t)=stats.tstat;
                    rec(t)=recovery;
                    t = t+1;
                end
            end
             if (i==12)
                for r=1:2
                    baseline_z = SUBzscore{i}{r}(lesion_intersect(GM_reshape==1)>thresh);
                    final_z = SUBzscore{i}{r+1}(lesion_intersect(GM_reshape==1)>thresh);
                    recovery = df(i,r+4)-df(i,r+3);
                    [h,p,ci,stats] = ttest2(baseline_z,final_z);
                    tst(t)=stats.tstat;
                    rec(t)=recovery;
                    t = t+1;
                end
            end
             if (i==20)
                for r=1
                    baseline_z = SUBzscore{i}{r}(lesion_intersect(GM_reshape==1)>thresh);
                    final_z = SUBzscore{i}{r+1}(lesion_intersect(GM_reshape==1)>thresh);
                    recovery = df(i,r+4)-df(i,r+3);
                    [h,p,ci,stats] = ttest2(baseline_z,final_z);
                    tst(t)=stats.tstat;
                    rec(t)=recovery;
                    t = t+1;
                end
             end
             if (i==22)
                for r=1:2
                    baseline_z = SUBzscore{i}{r+2}(lesion_intersect(GM_reshape==1)>thresh);
                    final_z = SUBzscore{i}{r+3}(lesion_intersect(GM_reshape==1)>thresh);
                    recovery = df(i,r+6)-df(i,r+5);
                    [h,p,ci,stats] = ttest2(baseline_z,final_z);
                    tst(t)=stats.tstat;
                    rec(t)=recovery;
                    t = t+1;
                end
             end
              if (i==23)
                for i=r:3
                    baseline_z = SUBzscore{i}{r+1}(lesion_intersect(GM_reshape==1)>thresh);
                    final_z = SUBzscore{i}{r+2}(lesion_intersect(GM_reshape==1)>thresh);
                    recovery = df(i,r+5)-df(i,r+4);
                    [h,p,ci,stats] = ttest2(baseline_z,final_z);
                    tst(t)=stats.tstat;
                    rec(t)=recovery;
                    t = t+1;
                end
            end
            for r=1:(nsess(i)-1)
                baseline_z = SUBzscore{i}{r}(lesion_intersect(GM_reshape==1)>thresh);
                final_z = SUBzscore{i}{r+1}(lesion_intersect(GM_reshape==1)>thresh);

                [h,p,ci,stats] = ttest2(baseline_z,final_z);
                tst(t)=stats.tstat;
                recovery = df(i,r+4)-df(i,r+3);
                rec(t)=recovery;
                t = t+1;
            end
        end
        
        plot(tst,rec, '.r', 'MarkerSize', 18) 
        hold on;
        b=polyfit(tst,rec,1);
        a=polyval(b,tst)
        plot(tst,a);

        [rho,p]=corr(tst', rec', 'Type', 'Pearson')
        pears_rho = rho
        pears_p=p
        ax =gca
        ax.FontSize=18
        txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        text(-90,-10,txt, 'fontsize', 18)
        xlabel('T-statistic: ICC between sequential sessions (Follow-up vs. Baseline)', 'fontsize', 18)
        ylabel('Sequential session change in Fugl-Meyer score (Follow-up - Baseline)', 'fontsize', 18)
        title({'Relationship between motor recovery & change in sequential baseline v. followup','ICC in cortical areas structurally connected to lesion'}, 'fontsize', 18)
        saveas(gcf, strcat(studydir, figuresdir, '/correlation_cortex_changeICCvschangeFuglMeyer_between-subsequent-FUs.png'))
    end
   % clear all;
    close all;

    % Figure 4 - change in ICC in areas structurally connected to the lesion specific regions (cerebellum, cortex, or brainstem)
     nvoxels=[];
     if ismember(figs, 4)
        for i=1:23
            if i==22
                continue
            end
            if i==23 
                continue
            end
            %load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
             
            % load motor recovery scores (final - initial)
            % column 8 = session 5
            % column 4 = session 1
            if i==6 % only 4 sessions (sessions 1, 2, 3, and 4)
                recovery = df(i,7)-df(i,4)
            elseif i==12 % only 3 sessions (session 1, 2, and 3)
                recovery = df(i,6)-df(i,4)
            elseif i==20 % only 2 sessions (session 1 and 2)
                recovery = df(i,5)-df(i,4)
            else
                recovery = df(i,8)-df(i,4)
            end
            disp(num2str(i));
            rec(i)=recovery
            k=1;
            mask = getMask(masks{k});
            mask = reshape(mask, [902629 1]);
            intersect = GM_reshape.*logical(mask);
            lesion_intersect = lesion_dc.*intersect;
   

            baseline_z = SUBzscore{i}{1}(lesion_intersect(GM_reshape==1) >thresh);  
            nvoxels(i,1) = size(baseline_z, 2);
            
            if i==6
                final_z = SUBzscore{i}{4}(lesion_intersect(GM_reshape==1)>thresh);
                histogram(baseline_z)
                hold on;
                histogram(final_z)
                title({'Z-scores at baseline (session 1) versus last follow up (session 4)', strcat('Num. voxels = ', num2str(nvoxels(i)))})
                xlabel('z-scored ICC')
                ylabel('count')
                saveas(gcf, strcat(studydir, figuresdir, 'distribution_S4_S1_SUB6.png'))
           clf;
            elseif i==12
                final_z = SUBzscore{i}{3}(lesion_intersect(GM_reshape==1)>thresh);
                histogram(baseline_z)
                hold on;
                histogram(final_z)
                title({'Z-scores at baseline (session 1) versus last follow up (session 3)', strcat('Num. voxels = ', num2str(nvoxels(i)))})
                xlabel('z-scored ICC')
                ylabel('count')
                saveas(gcf, strcat(studydir, figuresdir, 'distribution_S3_S1_SUB12.png'))
         clf;
            elseif i==20
                final_z = SUBzscore{i}{2}(lesion_intersect(GM_reshape==1)>thresh);
                histogram(baseline_z)
                hold on;
                histogram(final_z)
                title({'Z-scores at baseline (session 1) versus last follow up (session 2)', strcat('Num. voxels = ', num2str(nvoxels(i)))})
                xlabel('z-scored ICC')
                ylabel('count')
                saveas(gcf, strcat(studydir, figuresdir, 'distribution_S2_S1_SUB20.png'))
            clf;
            else % most subjects came for 5 sessions.
                final_z = SUBzscore{i}{5}(lesion_intersect(GM_reshape==1)>thresh);
                histogram(baseline_z)
                hold on;
                histogram(final_z)
                title({'Z-scores at baseline (session 1) versus last follow up (session 5)', strcat('Num. voxels = ', num2str(nvoxels(i)))})
                xlabel('z-scored ICC')
                ylabel('count')
                saveas(gcf, strcat(studydir, figuresdir, 'distribution_S5_S1_SUB', num2str(i), '.png'))
           clf;
            end
            [h,p,ci,stats]=ttest2(baseline_z,final_z,'Vartype','unequal');
            tst(i)=stats.tstat
            p(i)=p
   
        end
        plot(tst,rec, '.r', 'MarkerSize', 18) 
        hold on;
        b=polyfit(tst,rec,1);
        a=polyval(b,tst)
        plot(tst,a);

        [rho,p]=corr(tst', rec', 'Type', 'Pearson')
        pears_rho = rho
        pears_p=p
        [rho,p]=corr(tst', rec', 'Type', 'Spearman')
        spear_rho=rho
        spear_p=p
        xlabel('T-statistic: ICC Session 5 vs. ICC Sesssion 1', 'fontsize', 18)
        ylabel('Session 5 - Session 1 Fugl-Meyer score', 'fontsize', 18)
        ax = gca
        ax.FontSize = 18
        txt = ({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        %txt2 = ({strcat('Spearmans Rank Correlation = ', num2str(spear_rho)), strcat('p = ',num2str(round(spear_p,4)))})
        text(0,70,txt, 'fontsize', 18)
        title({'Relationship between motor recovery and change in session 1 vs session 5 ICC in ', masks{k}, ' & structurally connected to lesion'}, 'fontsize', 18)
        saveas(gcf, strcat(studydir, figuresdir, '/correlation_', masks{k}, '_changeICCvschangeFuglMeyer_baseline-vs-lastFU.png'))

end

% Figure 5 - correlation: change in ICC between session 1 and session2 & change in Fugl-Meyer session 1 and sesssion 2.
    if ismember(5,figs)
        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1000 700])

        tst_1=[]; %tstatistic (change in ICC)
        rec_1=[]; %recovery (change in fugl-meyer)
        tst_2=[];
        rec_2=[];
        tst_3=[];
        rec_3=[];
        tst_4=[];
        rec_4=[];

        disp(num2str(i))
        %load disconnectivity files
        lesion_dc=read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
        lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
        k=1;
        mask = getMask(masks{k});
        mask = reshape(mask, [902629 1]);
        intersect = GM_reshape.*logical(mask);
        lesion_intersect = lesion_dc.*intersect;
        
        % S2-S1
        t=1;
        for i=1:23
            if (i==23 | i==22)
                continue;
            end
            baseline_z = SUBzscore{i}{1}(lesion_intersect(GM_reshape==1)>thresh);
            final_z = SUBzscore{i}{2}(lesion_intersect(GM_reshape==1)>thresh);
            recovery = df(i,5)-df(i,4);
            [h,p,ci,stats] = ttest2(baseline_z,final_z);
            tst_1(t)=stats.tstat;
            rec_1(t)=recovery;
            t=t+1;
        end

        % S3-S2
        t=1;
        for i=1:23
            if (i==22 | i==20)
                continue
            end
            baseline_z = SUBzscore{i}{2}(lesion_intersect(GM_reshape==1)>thresh);
            final_z = SUBzscore{i}{3}(lesion_intersect(GM_reshape==1)>thresh);
            recovery = df(i,6)-df(i,5);
            [h,p,ci,stats] = ttest2(baseline_z,final_z);
            tst_2(t)=stats.tstat;
            rec_2(t)=recovery;
            t=t+1;
        end

        % S4-S3 
        for i=1:23 
            if (i==20 | i==12) 
                continue;
            end
            baseline_z = SUBzscore{i}{3}(lesion_intersect(GM_reshape==1)>thresh);
            final_z = SUBzscore{i}{4}(lesion_intersect(GM_reshape==1)>thresh);
            recovery = df(i,7)-df(i,6);
            [h,p,ci,stats] = ttest2(baseline_z,final_z);
            tst_3(t)=stats.tstat;
            rec_3(t)=recovery;
            t=t+1;
        end

        t=1;
        % S5-S4 
        for i=1:23 
            if ( i==6 | i==12 | i==20)
                continue;
            end
            baseline_z = SUBzscore{i}{4}(lesion_intersect(GM_reshape==1)>thresh);
            final_z = SUBzscore{i}{5}(lesion_intersect(GM_reshape==1)>thresh);
            recovery = df(i,8)-df(i,7);
            [h,p,ci,stats] = ttest2(baseline_z,final_z);
            tst_4(t)=stats.tstat;
            rec_4(t)=recovery;
            t=t+1;
        end

        subplot(2, 2, 1);    
        plot(tst_1,rec_1, '.r', 'MarkerSize', 18) 
        hold on;
        b=polyfit(tst_1,rec_1,1);
        a=polyval(b,tst_1)
        plot(tst_1,a);
        [rho,p]=corr(tst_1', rec_1', 'Type', 'Pearson')
        pears_rho = rho
        pears_p=p
        ax =gca
        ax.FontSize=18
        txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        text(-20,-20,txt, 'fontsize', 12)
            sgtitle('Session 1 vs Session 2')
        xlabel('T-statistic: ICC between Session 2 and Session 1', 'fontsize', 12)
        ylabel('Sequential session change in Fugl-Meyer score (Session 2 - Session 1)', 'fontsize', 12)
    
        subplot(2, 2, 2);
        plot(tst_2,rec_2, '.r', 'MarkerSize', 18) 
        hold on;
        b=polyfit(tst_2,rec_2,1);
        a=polyval(b,tst_2)
        plot(tst_2,a);
        [rho,p]=corr(tst_2', rec_2', 'Type', 'Pearson')
        pears_rho = rho
        pears_p=p
        ax =gca
        ax.FontSize=18
        txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        text(-20,-20,txt, 'fontsize', 12)
         sgtitle('Session 2 vs Session 3')
        xlabel('T-statistic: ICC between Session 3 and Session 2', 'fontsize', 12)
        ylabel('Sequential session change in Fugl-Meyer score (Session 3 - Session 2)', 'fontsize', 12)
    
        subplot(2, 2, 3);
        plot(tst_3,rec_3, '.r', 'MarkerSize', 18) 
        hold on;
        b=polyfit(tst_3,rec_3,1);
        a=polyval(b,tst_3)
        plot(tst_3,a);
        [rho,p]=corr(tst_3', rec_3', 'Type', 'Pearson')
        pears_rho = rho
        pears_p=p
        ax =gca
        ax.FontSize=18
        txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        text(-20,-20,txt, 'fontsize', 14)
        sgtitle('Session 3 vs Session 4')
        xlabel('T-statistic: ICC between Session 4 and Session 3', 'fontsize', 14)
        ylabel('Sequential session change in Fugl-Meyer score (Session 4 - Session 3)', 'fontsize', 14)

        subplot(2, 2, 4);
        plot(tst_4,rec_4, '.r', 'MarkerSize', 18) 
        hold on;
        b=polyfit(tst_4,rec_4,1);
        a=polyval(b,tst_4)
        plot(tst_4,a);
        [rho,p]=corr(tst_4', rec_4', 'Type', 'Pearson')
        pears_rho = rho
        pears_p=p
        ax =gca
        ax.FontSize=18
        txt=({strcat('Pearsons Correlation = ', num2str(pears_rho)), strcat('p = ',num2str(round(pears_p,4)))})
        text(-20,-20,txt, 'fontsize', 12)
         sgtitle('Session 5 vs Session 5')
        xlabel('T-statistic: ICC between Session 5 and Session 4', 'fontsize', 12)
        ylabel('Sequential session change in Fugl-Meyer score (Session 5 - Session 4)', 'fontsize', 12)
        saveas(gcf, strcat(studydir, figuresdir, '/correlation_cortex_changeICCvschangeFuglMeyer-separatedSessions.png'))
     end 

     %Figure 6 - plotting motor recovery information (Fugl Meyer)
     if ismember(6, figs)
        for i=1:23
            ran{i}=rand(1,3)
            % FUGL MEYER
            % session 1 = column 4, session 5 = column 8
            df(6, 8)=NaN; %no data
            df(12,7)=NaN;%no data
            df(12,8)=NaN;%no data
            df(20, 6:8)=NaN;%no data
            df(22, 4:5)=NaN;%no data
            df(23, 4)=NaN;%no data
        
            plot([7 14 30 90 180], df(i,4:8), '-o', 'color', ran{i}, 'LineWidth',1.5)
            title({'Fugl-Meyer scores of subject', num2str(i)})
            ylim([0 100])
            xlim([7 180])
            xticks([7 14 30 90 180])
            xticklabels({'7', '14', '30', '90', '180'})
            xlabel('Days post stroke', 'FontSize', 17)
            ylabel('Fugl-Meyer motor ability score', 'FontSize', 16)
            title(strcat('Subject ', num2str(i)))
            saveas(gcf, strcat(studydir, figuresdir,'fugl_meyer_scaled_space',num2str(i),'.png'));

            %GRIP STRENGTH
            % session 1=9,10 (left, right)
            % session 5=17,18 (left, right)
            df(6,17:18)=NaN; %no data
            df(12,15:18)=NaN;%no data
            df(20, 13:18)=NaN;%no data
            df(22, 9:12)=NaN;%no data
            df(23, 9:10)=NaN;%no data

            plot([7 14 30 90 180], df(i,9:2:17), '-o', 'color', ran{i}, 'LineWidth', 1.5)
            hold on;
            plot([7 14 30 90 180], df(i,10:2:18), '-o', 'color', ran{i+1}, 'LineWidth', 1.5)
            legend('Left', 'Right')
            ylim([0 100])
            xlim([7 180])
            xticks([7 14 30 90 180])
            xticklabels({'7', '14', '30', '90', '180'})
            xlabel('Days post stroke', 'FontSize', 17)
            ylabel('Fugl-Meyer motor ability score', 'FontSize', 16)
            if df(i, 2)==0
                lesionSide='Left'
            else
                lesionSide='Right'
            end
            if df(i,3)==0
                handedness='Left'
            else
                handedness='Right'
            end

            title({strcat('Subject ', num2str(i), ', stroke side =', lesionSide, ', handedness =', handedness)})
            saveas(gcf, strcat(studydir, figuresdir,'grip_strength',num2str(i),'.png'));

        end
    end


  %  elseif ismember(figs, 6)
  %  elseif ismember(figs, 7)
end
