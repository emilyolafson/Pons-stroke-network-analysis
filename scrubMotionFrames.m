function [outlierframes] = scrubMotionFrames(i, j, nscans, sub_direct)
    % Motion scrubbing - flag frames with high motion & 1 frame before + after
    subject = strcat('SUB', num2str(i));
    disp(strcat('flagging high-motion frames for:', subject))
    %load motion outlier files
    if i > 23
        outliers=load(strcat(sub_direct,subject,'/func/art_regression_outliers_aufunc1.mat')); %row=frames. one column for each flagged frame.
    else
        outliers=load(strcat(sub_direct,subject,'/func/S',num2str(j),'/art_regression_outliers_aufunc1.mat')); %row=frames. one column for each flagged frame.
    end
    
    outlierframes=outliers.R;

    motionframes=find(sum(outlierframes,2)); %returns frames which are high motion 
    motionframes=unique([motionframes; motionframes+1]); %add frame after high-motion frame (frame before is already flagged as an outlier by CONN).
    motionframes=motionframes(motionframes<=size(outlierframes,1)); %cut off frame 125 if added by above step.
    outlierframes=zeros(size(outlierframes,1), size(motionframes,1)); %outlierframes = empty matrix to be populated below.

    %remove first 5 frames of each scan.
    if nscans(j,i)==4
        motionframes=[motionframes;1;2;3;4;5;125;126;127;128;129;249;250;251;252;253;373;374;375;378;379];
    end
    if nscans(j,i)==3
        motionframes=[motionframes;1;2;3;4;5;125;126;127;128;129;249;250;251;252;253];
    end
    sort(unique(motionframes));

    %set cells to 1 if they are outliers (one outlier per row/column) (adds frame after ART-identified frame - the one before has already been set to 1, by ART)
    for l=1:size(motionframes,1)
        outlierframes(motionframes(l),l) = 1; 
    end
end
