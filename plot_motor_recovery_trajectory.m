%plot motor recovery of stroke patients
for i=1:23
    ran{i}=rand(1,3)
end

df = readmatrix(strcat(studydir, 'demog_strokepts2.csv'));

figure(1)
for i=1:23 %first row is column titles
    figure(i-1)
    plot([7 14 30 90 180], df(i,2:6), '-o', 'color', ran{i-1}, 'LineWidth',1.5)
    title({'Fugl-Meyer scores of subject', num2str(i)})
    ylim([0 100])
    xlim([7 180])
    xticks([7 14 30 90 180])
    xticklabels({'7', '14', '30', '90', '180'})
    xlabel('Days post stroke', 'FontSize', 17)
    ylabel('Fugl-Meyer motor ability score', 'FontSize', 16)
    title(strcat('Subject ', num2str(i-1)))
    saveas(gcf, strcat(studydir, figuresdir,'fugl_meyer_scaled_space',num2str(i-1),'.png'));
end

a=get(gca, 'XTickLabel')
%legend({'SUB1', 'SUB2', 'SUB3', 'SUB4', 'SUB5', 'SUB6', 'SUB7', 'SUB8', 'SUB9', 'SUB10', 'SUB11'})

set(gca,'XTickLabel', a, 'fontsize', 18)
set(gcf, 'Position',[0 0 1200 1200])

figure(2)
gripstrength=readcell('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_gripstrength.csv')
for i=2:12
    if i==7
        continue
    end
    if length(setdiff(gripstrength(i,3), {'L'}))==0 %left handed, right cortex dominant.
        plot(cell2mat(gripstrength(i, [4 6 8 10 12])), '-o', 'color', ran{i-1}, 'LineWidth', 1.5)
        hold on;
    end
    
    if length(setdiff(gripstrength(i,3), {'R'}))==0
        plot(cell2mat(gripstrength(i, [5 7 9 11 13])), '-o', 'color', ran{i-1}, 'LineWidth', 1.5)
        hold on;
    end
end
xticks([7 14 30 90 180])

ylim([0 100])
xlim([7 180])
xticks([7 14 30 90 180])
xticklabels({'7', '14', '30', '90', '180'})
xlabel('Days post stroke', 'FontSize', 17)
ylabel('Fugl-Meyer motor ability score', 'FontSize', 16)
a=get(gca, 'XTickLabel')
set(gca,'XTickLabel', a, 'fontsize', 18)




figure(3)%barthel index
df =readmatrix('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_barthel_index.csv')
df(7,6) = NaN
for i=2:12 %first row is column titles
    figure(i-1)
    plot([7 14 30 90 180], df(i,2:6), '-o', 'color', ran{i-1}, 'LineWidth',1.5)
    title('Barthel Index scores of 11 stroke patients (higher = more recovery)')
    ylim([0 100])
    xlim([7 180])
    xticks([7 14 30 90 180])
    xticklabels({'7', '14', '30', '90', '180'})
    xlabel('Days post stroke', 'FontSize', 17)
    ylabel('Barthel Index', 'FontSize', 16)
    hold on;
end  
hold off;

a=get(gca, 'XTickLabel')
set(gca,'XTickLabel', a, 'fontsize', 18)
set(gcf, 'Position',[0 0 1200 1200])


figure(4)%NIHSS
df =readmatrix('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_NIHSS.csv')
df(7,6) = NaN
for i=2:12 %first row is column titles
    plot([7 14 30 90 180], df(i,2:6), '-o', 'color', ran{i-1}, 'LineWidth',1.5)
    title('NIHSS scores of 11 stroke patients (lower = more recovery)')
    ylim([0 6])
    xlim([7 180])
    xticks([7 14 30 90 180])
    xticklabels({'7', '14', '30', '90', '180'})
    xlabel('Days post stroke', 'FontSize', 17)
    ylabel('NIHSS', 'FontSize', 16)
    hold on;
end  
hold off;
a=get(gca, 'XTickLabel')
set(gca,'XTickLabel', a, 'fontsize', 18)
set(gcf, 'Position',[0 0 1200 1200])


figure(5)%modified rankin scale
df =readmatrix('/mnt/shared_data3/emo4002/stroke_pts/Demographics_11pts_mRS.csv')
df(7,6) = NaN
for i=2:12 %first row is column titles
    plot([7 14 30 90 180], df(i,2:6), '-o', 'color', ran{i-1}, 'LineWidth',1.5)
    title('Modified rankin scale scores of 11 stroke patients (lower = more recovery)')
    ylim([0 6])
    xlim([7 180])
    xticks([7 14 30 90 180])
    xticklabels({'7', '14', '30', '90', '180'})
    xlabel('Days post stroke', 'FontSize', 17)
    ylabel('modified rankin scale', 'FontSize', 16)
    hold on;
end  
hold off;
a=get(gca, 'XTickLabel')
set(gca,'XTickLabel', a, 'fontsize', 18)
set(gcf, 'Position',[0 0 1200 1200])

