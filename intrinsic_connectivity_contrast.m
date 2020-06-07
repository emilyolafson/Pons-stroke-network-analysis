function [ICC] = intrinsic_connectivity_contrast(timeSeries)
% [ICC] = INTRINSIC_CONNECTIVITY_CONTRAST(TIMESERIES)

% from here: 
% https://www.mathworks.com/matlabcentral/fileexchange/68248-intrinsic_connectivity_contrast

%   The intrinsic connectivity contrast (ICC) is a measure of degree 
%   centrality for weighted networks. The ICC for a given voxel i is the 
%   average of the squared correlation values of the voxel i time series 
%   with that of every other voxel (i.e., ICC = mean(r.^2)). This function 
%   enables rapid computation of the ICC on extremely large resting-state 
%   fMRI datasets, including cases in which explict calculation of the full 
%   correlation matrix is prohibitive, whether due to RAM or computation 
%   time limitations. 
% 
% Inputs:
%   timeSeries,   size [times x voxels] matrix
% 
% Outputs:
%   ICC,          a vector [1 x voxels] of Intrinsic Connectivity Contrast
%                 values
% 
% *Intrinsic Connectivity Contrast References:
% 
% Martuzzi, R., Ramani, R., Qiu, M., Shen, X., Papademetris, X., & 
%   Constable, R. T. (2011). A whole-brain voxel based measure of intrinsic 
%   connectivity contrast reveals local changes in tissue connectivity with 
%   anesthetic without a priori assumptions on thresholds or regions of 
%   interest. Neuroimage, 58(4), 1044-1050.
% 
% Whitfield-Gabrieli, S., & Nieto-Castanon, A. (2012). Conn: a functional 
%   connectivity toolbox for correlated and anticorrelated brain networks. 
%   Brain connectivity, 2(3), 125-141.
% 
% For implementation details, see Whitfield-Gabrieli & Nieto-Castanon (2012)
%   Appendix A.2: Computing Voxel-Wise Linear Functional Connectivity 
%   Measures Singular Value Decomposition
% 
% Author:  Elliot Layden, The University of Chicago, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nTimes = size(timeSeries,1); nVox = size(timeSeries,2);
    X = timeSeries; 
    X = bsxfun(@minus,X,nansum(X,1)./nTimes); % mean center
    X = X.*repmat(sqrt(1./max(eps,nansum(abs(X).^2,1))),[nTimes,1]);  % L2 norm
    [U,S,~] = svd(cov(X')); % time x time cov matrix
    S = diag(S); b = U'*X; 
    ICC = sum(repmat(S,1,nVox).*b.^2);
end
