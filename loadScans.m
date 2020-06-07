function cat_func = loadScans(i, j, nscans, sub_direct)
  % loadScans: loads denoised fMRI scans from a single session concatenated together across time.
    cat_func=[];
   
    subject = strcat('SUB', num2str(i));
    disp(strcat('Loading scans for subejct:', num2str(i), ', session:', num2str(j)));

    for k=1:nscans(j,i) %loop over scans within sessions
      if i > 23 %controls don't have session directories
          func=read_avw(strcat(sub_direct,subject,'/func/denoise_swaufunc',num2str(k),'.nii')); %91x109x91
        else
          func=read_avw(strcat(sub_direct,subject,'/func/S',num2str(j),'/denoise_swaufunc',num2str(k),'.nii')); %91x109x91
        end
        func=reshape(func,[], size(func,4)); %flattened 2D matrix that is <voxels>
        if nscans(j,i)==1
            cat_func=func;
        else
          cat_func=horzcat(cat_func, func); % concatenate scans together in time. should have 124*#scans columns
        end 
    end
end
