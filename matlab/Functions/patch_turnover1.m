function [patchID,land]=patch_turnover1(matsize,patchID,n_patch,pturnover)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function created to simulate patch turnover in a patchy landscape
% (patchsize=1)
%
% input is:
% matsize (1D size of matrix)
% patchID( matrix with patches with unique identifyers)
% n_patch  ( vector with unique patch IDs)
% pturnover (fraction of patches that is relocated)
%
% output is:
% land ( matrix identifying habitat and non habitat grid cells )
% patchID( matrix with patches with unique identifyers)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


id_pt=unique(n_patch(rand(length(n_patch),1) < pturnover));
patchID(ismember(patchID,id_pt))=0;   
for i = 1:length(id_pt)
    % randomly select midpoint of new patch
    x = randi(matsize); y = randi(matsize);
    
    % accept or reject based on distance
    while patchID(y,x)>0
        x = randi(matsize); y = randi(matsize);
    end
    patchID(y,x)=id_pt(i);
end

land = 1*(patchID>0);