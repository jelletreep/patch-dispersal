function [patchID,n_patch,land]=land_gen(matsize,numpatch)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function created to generate landscape matrix with circular patches
%
% input is:
% matsize (1D size of matrix)
% numpatch (total number of patches in domain)
%
% output is:
% patchID (landscape matrix with patches with unique IDs)
% n_patch  ( vector with unique patch IDs)
% land ( matrix identifying habitat and non habitat grid cells )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



convland = zeros(matsize);
for i = 1:numpatch
    % randomly select midpoint of new patch
    x = randi(matsize); y = randi(matsize);
    
    % accept or reject based on distance
    while convland(y,x)>0
        x = randi(matsize); y = randi(matsize);
    end
    
    convland(y,x)=i;
    
end
patchID = convland;
n_patch = unique(patchID); n_patch=n_patch(2:end);
land = 1*(patchID>0);
