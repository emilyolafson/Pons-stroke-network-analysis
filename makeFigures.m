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

    % Figure 1 - boxplots of voxels connected vs unconnected to lesion

    nsess=[5;5;5;5;5;4;5;5;5;5;5;3;5;5;5;5;5;5;5;2;5;5;5];

    thresh = 0.1;
    if ismember(figs, 1)
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
        
    elseif ismember(figs, 2) 
    % Figure 2 - 
    elseif ismember(figs, 3)
    elseif ismember(figs, 4)
    elseif ismember(figs, 5)
    elseif ismember(figs, 6)
    elseif ismember(figs, 7)
    end

en

