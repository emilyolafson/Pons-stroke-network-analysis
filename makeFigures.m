function [] = makeFigures(nsess, studydir, figuresdir, disconnectivitydir. figs=1:5)
% makeFigures: generate figures for Pontine stroke network analysis project. 
% INPUT: 
%     nsess: number of sessions per subject (columns)
%     studydir: folder containing subject folders and folder for outputs.
%     figures: name of folder where figures will be saved.
%     disconnectivitydir: name of folder where disconnectivity maps are.
%     figs: list of figures to be generated.
%
    % load z-scores
    load(strcat(studydir, resultsdir, 'zICC_23subjects.mat'),'SUBzscore'))
    
    % load gray matter masks
    GM = read_avw('/home/emo4002/colossus_shared3/c1referenceT1.nii'); 
    GM_reshape = reshape(GM, [902629 1]);
    GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
    GM_reshape(GM_reshape <= 0.25) = 0;

    % set analysis variables
    thresh = 0.1; % value above which voxels are considered to be 'connected' to the lesion

    % Figure 1 - boxplots of voxels connected vs unconnected to lesion
    if ismember(figs, 1)
        k=0;
        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1400 1400])
        p = 1;
        for i = 1:5
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(5,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig1)

        fig2 = figure(2)
        set(fig1, 'Position', [0 0 1400 280])
        p = 1;
        for i = 6
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(1,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig2)

        fig3 = figure(3)
        set(fig3, 'Position', [0 0 1400 1400])
        p = 1;
        for i = 7:11
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(5,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig3)

        fig4 = figure(4)
        set(fig4, 'Position', [0 0 1400 280])
        p = 1;
        for i = 12
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(1,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig4)

        fig5 = figure(5)
        set(fig5, 'Position', [0 0 1400 1400])
        p = 1;
        for i = 13:17
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(5,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig5)

        fig5 = figure(5)
        set(fig5, 'Position', [0 0 1400 560])
        p = 1;
        for i = 18:19
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(2,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig6)

        fig6 = figure(6)
        set(fig6, 'Position', [0 0 1400 280])
        p = 1;
        for i = 20
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(1,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig6)

        fig7 = figure(7)
        set(fig7, 'Position', [0 0 1400 840])
        p = 1;
        for i = 21:23
            subject = strcat('SUB',num2str(i));
            % load disconnectivity files
            lesion_dc = read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc = reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
            for j = 1:nsess(i)
                pos = SUBzscore{i,j}(lesion_dc(GM_reshape)>thresh);
                sizepos = size(pos,2);
                zer = SUBzscore{i,j}(lesion_dc(GM_reshape)<thresh);
                sizezer = size(zer,2);
                grp = [zeros(1,size(pos,2)), ones(1,size(zer,2))];
                fun = [pos zer];
                hold on;
                subplot(3,5, p);
                boxplot(fun, grp);
                ylim([-6 6])
                title({strcat('SUB', num2str(i), ' session ', num2str(j)) , strcat("n=", num2str(sizepos), ", n=", num2str(sizezer))})
                xlabel('Connectivity to lesion area')
                ylabel('z-score')
                xticklabels({'Connected', 'Not connected'})
                p=p+1;
            end
        end
        k=k+1;
        exportgraphics(gcf, strcat(studydir, figuresdir, '/boxplots_connect-vs-disconnect_', num2str(k), '.pdf', 'ContentType', 'vector');
        close(fig7)
    end
    clear all;
    close all;


    % Figure 2 - change in motor scores vs change in ICC from the last scan to the first scan 
    if  ismember(figs, 2) 
        df = readmatrix(studydir, 'Demographics_ponsstroke_23subjs.csv');
        tst=[]; %tstatistic (change in ICC)
        rec=[]; %recovery (change in fugl-meyer)

        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1400 1400]);

        for i=1:23
            subject = strcat('SUB',num2str(i));
            
            %load gray matter masks
            GM = read_avw(str_cat(studydir, 'c1referenceT1.nii'));
            GM_reshape = reshape(GM, [1 902629])
            GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
            GM_reshape(GM_reshape <=0.25) = 0;

            %load disconnectivity files
            lesion_dc=read_avw(strcat(studydir, disconnectivitydir, SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
                
            baseline_z = SUBzscore{i,1}(lesion_dc(GM_reshape)>thresh);

            if i==6
                final_z = SUBzscore{i,4}(lesion_dc(GM_reshape)>thresh);
            elseif i==12
                final_z = SUBzscore{i,3}(lesion_dc(GM_reshape)>thresh);
            elseif i==20
                final_z = SUBzscore{i,2}(lesion_dc(GM_reshape)>thresh);
            else % most subjects came for 5 sessions.
                final_z = SUBzscore{i,5}(lesion_dc(GM_reshape)>thresh);
            end
            
            % calculate the difference in z-score normalized ICC values at baseline and followup.
            % metric is a t-statistic using a paired t-test with unequal variances.
            % positive t-stat means distribution 
            [h,p,ci,stats]=ttest2(baseline_z,follow_up,'Vartype','unequal')
            tst(i)=stats.tstat
            p(i)=p

            % load motor recovery scores (final - initial)
            % column 6 = sesssion 1
            % column 10 = sesssion 5
            if i==6 % only 4 sessions (sessions 1, 2, 3, and 4)
                recovery = df(i+3,9)-df(i+3,6)
            elseif i==12 % only 3 sessions (session 1, 2, and 3)
                recovery = df(i+3,8)-df(i+3,6)
            elseif i==20 % only 2 sessions (session 1 and 2)
                recovery = df(i+3,7)-df(i+3,6)
            else
                recovery = df(i+3,10)-df(i+3,6)
            end
            rec(i)=recovery
            savea(gcf, )
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
        exportgraphics(gcf, strcat(studydir, figuresdir, '/correlation_allGM_changeICCvschangeFuglMeyer_baseline-vs-lastFU.pdf', 'ContentType', 'vector'))
    end
    clear all;
    close all;

    % Figure 3 - Boxplots of each inter-session change in ICC
    if ismember(figs, 3)
        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1400 1400])

        tst=[]; %tstatistic (change in ICC)
        rec=[]; %recovery (change in fugl-meyer)

        %load gray matter masks
        GM = read_avw(str_cat(studydir, 'c1referenceT1.nii'));
        GM_reshape = reshape(GM, [1 902629])
        GM_reshape(GM_reshape > 0.25) = 1; %threshold GM mask
        GM_reshape(GM_reshape <=0.25) = 0; lesion_dc{w} = read_avw(strcat('/home/emo4002/colossus_shared3/stroke_pts/',subject,'/lesionTract/HCP_MGH_32fold/wc_',subject,'_lesionTract.nii.gz'));
        df = readmatrix(studydir, 'Demographics_ponsstroke_23subjs.csv');

        t=1;
        for i=1:23
            %load disconnectivity files
            lesion_dc=read_avw(strcat(studydir, disconnectivitydir, SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
                
            for r=1:(nsess(i)-1)
                baseline_z = SUBzscore{i,r}(lesion_dc(GM_reshape)>thresh);
                final_z = SUBzscore{i,r+1}(lesion_dc(GM_reshape)>thresh);

                [h,p,ci,stats] = ttest2(baseline_z,final_z);
                tst(t)=stats.tstat
                recovery = df(w+1,r+2)-df(w+1,r+1)
                rec(t)=recovery
                t = t+1;
            end
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
        xlabel('T-statistic: ICC between sequential sessions (Follow-up vs. Baseline)', 'fontsize', 18)
        ylabel('Sequential session change in Fugl-Meyer score (Follow-up - Baseline)', 'fontsize', 18)
        title({'Relationship between motor recovery & change in sequential baseline v. followup','ICC in cortical areas structurally connected to lesion'}, 'fontsize', 18)
        exportgraphics(gcf, strcat(studydir, figuresdir, '/correlation_allGM_changeICCvschangeFuglMeyer_between-subsequent-FUs.pdf', 'ContentType', 'vector'))
    end
    clear all;
    close all;

    % Figure 4 - change in ICC in areas structurally connected to the lesion specific regions (cerebellum, cortex, or brainstem)
    if ismember(figs, 4)
    end


    % Figure 4 - change in ICC in areas structurally connected to the lesion specific regions (cerebellum, cortex, or brainstem)
    elseif ismember(figs, 5)
    elseif ismember(figs, 6)
    elseif ismember(figs, 7)
    end

en

