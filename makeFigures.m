function [] = makeFigures(nsess, studydir, figuresdir, resultsdir, disconnectivitydir, figs)
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
    df = readmatrix(strcat(studydir, 'Demographics_PonsStroke_23subs.csv'));

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

    % Figure 3 - Boxplots of each inter-session change in ICC
    if ismember(3, figs)
        fig1 = figure(1)
        set(fig1, 'Position', [0 0 1000 1000])

        tst=[]; %tstatistic (change in ICC)
        rec=[]; %recovery (change in fugl-meyer)

        t=1;
        for i=1:23
            %load disconnectivity files
            lesion_dc=read_avw(strcat(studydir, disconnectivitydir, 'SUB',num2str(i), '_voxeldisconnect_2mm.nii.gz'));
            lesion_dc=reshape(lesion_dc,[902629 1]); %flattened 1D matrix that is <voxels>
                
            for r=1:(nsess(i)-1)
                baseline_z = SUBzscore{i}{r}(lesion_dc(GM_reshape)>thresh);
                final_z = SUBzscore{i}{r+1}(lesion_dc(GM_reshape)>thresh);

                [h,p,ci,stats] = ttest2(baseline_z,final_z);
                tst(t)=stats.tstat
                recovery = df(i+3,r+7)-df(i+3,r+6)
                rec(t)=recovery
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
        text(-90,-60,txt, 'fontsize', 18)
        xlabel('T-statistic: ICC between sequential sessions (Follow-up vs. Baseline)', 'fontsize', 18)
        ylabel('Sequential session change in Fugl-Meyer score (Follow-up - Baseline)', 'fontsize', 18)
        title({'Relationship between motor recovery & change in sequential baseline v. followup','ICC in cortical areas structurally connected to lesion'}, 'fontsize', 18)
        saveas(gcf, strcat(studydir, figuresdir, '/correlation_allGM_changeICCvschangeFuglMeyer_between-subsequent-FUs.png'))
    end
   % clear all;
    close all;

    % Figure 4 - change in ICC in areas structurally connected to the lesion specific regions (cerebellum, cortex, or brainstem)
    if ismember(figs, 4)
    end


    % Figure 4 - change in ICC in areas structurally connected to the lesion specific regions (cerebellum, cortex, or brainstem)
  %  elseif ismember(figs, 5)
  %  elseif ismember(figs, 6)
  %  elseif ismember(figs, 7)
end
