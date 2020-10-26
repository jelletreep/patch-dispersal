function [n_patch,land,convland]=land_gen2(matsize,patchsize,numpatch)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function created to generate landscape matrix with circular patches
%
% input is:
% matsize (1D size of matrix)
% patchsize (diameter of patch in number of gridcells)
% numpatch (total number of patches in domain)
%
% output is:
% n_patch  ( vector with unique patch IDs)
% land ( matrix identifying habitat and non habitat grid cells )
% convland ( matrix with extended boundaries and patches with unique
% identifyers)
%
% Domain has periodic boundaries. Patches can cross border of periodic
% boundaries.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if patchsize == 2
    convland = zeros(matsize+3*patchsize);
    
    for i = 1:numpatch
        
        % randomly select midpoint of new patch
        x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
        % accept or reject based on distance
        while sum(sum(convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))))>0
            x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
        end
        
        convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=i;
        
        % account for periodic boundaries
        if x <= 2.5*patchsize
            convland((y-patchsize/2):(y+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=i;
        elseif x > matsize+0.5*patchsize
            convland((y-patchsize/2):(y+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=i;
        end
        if y <= 2.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=i;
        elseif y > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=i;
        end
        if y <= 2.5*patchsize && x <= 2.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=i;
        elseif y <= 2.5*patchsize && x > matsize+0.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=i;
        elseif x <= 2.5*patchsize && y > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=i;
        elseif y > matsize+0.5*patchsize && x > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=i;
        end
        
    end
    patchID=convland((1.5*patchsize+1):(1.5*patchsize+matsize),(1.5*patchsize+1):(1.5*patchsize+matsize));
    
elseif patchsize > 2
    
    i=1;
    [cc,rr] = meshgrid(1:(matsize+3*patchsize));
    convland = zeros(matsize+3*patchsize);
    
    while i <= numpatch
        
        it = 0;
        ext = 0;
        
        % randomly select midpoint of new patch
        x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
        % accept or reject based on distance
        
        while sum(sum(convland((y-patchsize/2):(y+patchsize/2-1),...
                (x-patchsize/2):(x+patchsize/2-1)).*...
                (sqrt((rr((y-patchsize/2):(y+patchsize/2-1),...
                (x-patchsize/2):(x+patchsize/2-1))+0.5-y).^2+...
                (cc((y-patchsize/2):(y+patchsize/2-1),...
                (x-patchsize/2):(x+patchsize/2-1))+0.5-x).^2)<=0.5*patchsize)))>0 && ext == 0;
            x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
            
            it = it+1;
            
            %% exit loop and restart when iterations exceed 50000
            if it > 50000
                convland = zeros(matsize+3*patchsize);   % possible without or statement???? otherwise replace with for loop
                i=1;
                ext = 1;
            end
            
        end
        
        C = sqrt((rr+0.5-y).^2+(cc+0.5-x).^2);
        convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=...
            convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))+...
            (C((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))<0.5*patchsize)*i;
        
        % account for periodic boundaries
        if x <= 2.5*patchsize
            convland((y-patchsize/2):(y+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=0;
            convland((y-patchsize/2):(y+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=...
                convland((y-patchsize/2):(y+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        elseif x > matsize+0.5*patchsize
            convland((y-patchsize/2):(y+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=0;
            convland((y-patchsize/2):(y+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        end
        if y <= 2.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=0;
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=...
                convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        elseif y > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=0;
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=...
                convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        end
        if y <= 2.5*patchsize && x <= 2.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=0;
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=...
                convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        elseif y <= 2.5*patchsize && x > matsize+0.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=0;
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=...
                convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        elseif x <= 2.5*patchsize && y > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=0;
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=...
                convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        elseif y > matsize+0.5*patchsize && x > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=0;
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=...
                convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))+...
                convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1));
        end
        i=i+1;
        
    end
    patchID=convland((1.5*patchsize+1):(1.5*patchsize+matsize),(1.5*patchsize+1):(1.5*patchsize+matsize));
end
n_patch = unique(patchID); n_patch=n_patch(2:end);
land = 1*(patchID>0);