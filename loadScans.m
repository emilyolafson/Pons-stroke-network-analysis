function cat_func = loadScans(i, j, nscans)
% Save scans from a single session concatenated together across time.
    cat_func=[];
   
    subject = strcat('SUB', num2str(i));
    for k=1:nscans(j,i) %loop over scans within sessions
        disp(strcat('loading functional scans for:', subject, ', session:', num2str(j), ',scan:', num2str(k)))
        func=read_avw(strcat('/home/emo4002/colossus_shared3/pons_sfmodelling/stroke_pts/',subject,'/func/S',num2str(j),'/denoise_swaufunc',num2str(k),'.nii')); %91x109x91
        func=reshape(func,[], size(func,4)); %flattened 2D matrix that is <voxels>
        if nscans(j,i)==1
            cat_func=func;
        else
          cat_func=horzcat(cat_func, func); % concatenate scans together in time. should have 124*#scans columns
        end 
    end
end
