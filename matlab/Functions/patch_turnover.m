function [land,convland]=patch_turnover(pturnover,matsize,patchsize,convland,n_patch)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function created to simulate patch turnover in a patchy landscape with
% circular patches
%
% input is:
% pturnover (fraction of patches that is relocated)
% matsize (1D size of matrix)
% patchsize (diameter of patch in number of gridcells)
% convland ( matrix with extended boundaries and patches with unique identifyers)
% n_patch  ( vector with unique patch IDs)
%
% output is:
% land ( matrix identifying habitat and non habitat grid cells )
% convland ( matrix with extended boundaries and patches with unique identifyers)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

id_pt=unique(n_patch(rand(length(n_patch),1) < pturnover));

%% remove old patches
convland(ismember(convland,id_pt))=0;   % possible without or statement???? otherwise replace with for loop

%% relocate patch
if patchsize == 2
    i=1;
    while i <= length(id_pt)
        it = 0;
        ext = 0;
        % randomly select midpoint of new patch
        x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
        
        %       % accept or reject based on distance
        while sum(sum(convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))))>0 && ext == 0
            x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
            it = it+1;
            
            %% exit loop and restart when iterations exceed 50000
            if it > 50000
                convland(ismember(convland,id_pt))=0;   % possible without or statement???? otherwise replace with for loop
                i=1;
                ext = 1;
            end
        end
        
        convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=id_pt(i);
        
        % account for periodic boundaries
        if x <= 2.5*patchsize
            convland((y-patchsize/2):(y+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=id_pt(i);
        elseif x > matsize+0.5*patchsize
            convland((y-patchsize/2):(y+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=id_pt(i);
        end
        if y <= 2.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=id_pt(i);
        elseif y > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=id_pt(i);
        end
        if y <= 2.5*patchsize && x <= 2.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=id_pt(i);
        elseif y <= 2.5*patchsize && x > matsize+0.5*patchsize
            convland((y+matsize-patchsize/2):(y+matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=id_pt(i);
        elseif x <= 2.5*patchsize && y > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x+matsize-patchsize/2):(x+matsize+patchsize/2-1))=id_pt(i);
        elseif y > matsize+0.5*patchsize && x > matsize+0.5*patchsize
            convland((y-matsize-patchsize/2):(y-matsize+patchsize/2-1),(x-matsize-patchsize/2):(x-matsize+patchsize/2-1))=id_pt(i);
        end
        i=i+1;
    end
    patchID=convland((1.5*patchsize+1):(1.5*patchsize+matsize),(1.5*patchsize+1):(1.5*patchsize+matsize));
    
elseif patchsize > 2
    i=1;
    [cc,rr] = meshgrid(1:(matsize+3*patchsize));
    
    while i <= length(id_pt)
        it = 0;
        ext = 0;
        % randomly select midpoint of new patch
        x = randi(matsize)+1.5*patchsize; y = randi(matsize)+1.5*patchsize;
        
        % reject new midpoint of patch when there is overlap with different patch
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
                convland(ismember(convland,id_pt))=0;   % possible without or statement???? otherwise replace with for loop
                i=1;
                ext = 1;
            end
        end
        
        C = sqrt((rr+0.5-y).^2+(cc+0.5-x).^2);
        convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))=...
            convland((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))+...
            (C((y-patchsize/2):(y+patchsize/2-1),(x-patchsize/2):(x+patchsize/2-1))<0.5*patchsize)*id_pt(i);
        
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

land = 1*(patchID>0);


