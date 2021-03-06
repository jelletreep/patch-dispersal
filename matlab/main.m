%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code
%
% Supplement to the manuscript:
% 'Seed dispersal as a search strategy: dynamic and fragmented landscapes 
% select for multi-scale movement strategies in plants'
% Authors: Jelle Treep, Monique de Jager, Frederic Bartumeus,
% Merel B. Soons
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
addpath('./Functions')

%% User Input Variables
% Landscape variables
matsize = 128;                          % domain size; must be 2^n
patchsize = 16;                         % patchsize is the diameter of a circular patch measured in number of grid cells. possible values: 2^n, n=0:7
ipatchdist = 100;                        % inter patch distance (Mean free path - patch radius). possible values within range: [patchsize - 1000]
pturnover = 0.1;                        % patch disappears after dispersal and is replaced by a new patch at a different location. possible values [0-1]
mu = [2.5 , 3.5];                         % shape parameter of pareto distribution shaped dispersal kernel; 1 < mu <= 3
ntime = 200;                            % total number of timesteps

if ipatchdist < patchsize
    error('Error: "ipatchdist" must be greater than or equal to "patchsize"')
end

%% Landscape parameters
mfp = ipatchdist + patchsize/2;                         % mean free path
numpatch = round(matsize^2/(mfp*patchsize));            % number of patches
hcov = (numpatch*(pi*((patchsize/2)^2)))/(matsize^2);   % habitat cover

%% dispersal parameters
pseed = [100, 100];                                     % number of seeds of species 1 and 2
Xmax = 10000;                                           % truncation distance
xm = 0.5;                                               % minimum distance pareto distribution

%% time control parameters
imtime = 1;                                             % number of timesteps between graph updates
im_i = 0;                                               % image iterator
idx=1;

%% input parameters for landscape matrices
dx = 1; dy = 1;                                             % cellsize
midpoint = matsize/2+1;                                     % define midpoint of landscape domain
x = linspace(-(matsize*dx)/2,(matsize*dx)/2-dx,matsize);    % x axis of domain
y = linspace(-(matsize*dy)/2,(matsize*dy)/2-dy,matsize);    % y axis of domain
[xx,yy] = meshgrid(x,y);                                    % domain coordinates (relative to midpoint)
distmat = sqrt(xx.^2+yy.^2);                                % distance matrix (distance to midpoint)

%% create initial landscape
if patchsize > 1 % then an extra matrix (convland) is created to account for periodic boundaries
    [n_patch,land,convland]=land_gen2(matsize,patchsize,numpatch);
else
    [patchID,n_patch,land]=land_gen(matsize,numpatch);
end
%% population variables
SP = mat_gen(matsize,2);                                % random placement in matrix with equal probabilities
SP(land==0) = 0;                                        % remove individuals located in 'no habitat'
cumsp1 = zeros(1,ntime); cumsp2 = zeros(1,ntime);       % create variable to keep track of population history
ef = 1;                                                 % extinction fraction. 1 means all species die after dispersal (semelparous)

%% Visualize initial population distributions
imagesc(SP)
colormap([1 1 1
    0.8500 0.3250 0.0980
    0 0.4470 0.7410]);
cbh=colorbar('v');
caxis([-0.5 2.5])
set(cbh,'YTick',0:2)
set(cbh,'YTickLabel',{'No Habitat','Species 1','Species 2'})

%% Initialize dispersal kernels of species
K1 = 1/(2*pi)*((-mu(1)+2.0001)/(Xmax.^(-mu(1)+2.0001)-xm^(-mu(1)+2.0001))).*distmat.^(-mu(1));
K1(midpoint,midpoint) = 0; K1 = K1/sum(sum(K1));        % set midpoint to zero and normalize
K2 = 1/(2*pi)*((-mu(2)+2.0001)/(Xmax.^(-mu(2)+2.0001)-xm^(-mu(2)+2.0001))).*distmat.^(-mu(2));
K2(midpoint,midpoint) = 0; K2 = K2/sum(sum(K2));        % set midpoint to zero and normalize

%% Fast fourier transform kernels
fK1 = fft2(K1);
fK2 = fft2(K2);

%% create figure for dynamic visualization
fig=figure('Position',[100,100,1200,400]);

tic

%% dynamic loop
i=1;                % timestep
convergence = 0;    % alternative exit of while loop if population have stabilized

while i < ntime && convergence == 0
    
    %% Species states
    SP1 = ones(matsize).*(SP == 1);                         % Distribution of species 1
    SP2 = ones(matsize).*(SP == 2);                         % Distribution of species 2
    cumsp1(i) = sum(sum(SP1));  cumsp2(i) = sum(sum(SP2));  % save population states
    
    
    %% Disperse seeds
    fSP1 = fft2(SP1);       fSP2 = fft2(SP2);               % fast fourier transform of population distributions
    fp1 = fSP1.*fK1;        fp2 = fSP2.*fK2;                % Product of transforms
    p1 = pseed(1).*real( fftshift( ifft2( fp1 ) ) );        % Evaluate convolution (Population dispersal kernel)
    p2 = pseed(2).*real( fftshift( ifft2( fp2 ) ) );        % Evaluate convolution (Population dispersal kernel)
    
    p_norm = p1/(100*cumsp1(i));                            % normalize kernel
    disc_seed=mnrnd(round(100*cumsp1(i)),p_norm(:));        % multinomial sampling of seeds
    d1 = reshape(disc_seed,[matsize,matsize]);              % create seed distribution
    
    p_norm = p2/(100*cumsp2(i));                            % normalize kernel
    disc_seed=mnrnd(round(100*cumsp2(i)),p_norm(:));        % multinomial sampling of seeds
    d2 = reshape(disc_seed,[matsize,matsize]);              % create seed distribution
    
    %% Species extinction
    ext_mat = rand(matsize) > ef;                           % determine which individuals die
    SP = SP.*ext_mat;                                       % Update population matrix
    
    %% Update landscape  (patch turnover; relocates patches randomly)
    if pturnover == 1 && patchsize > 1
        [n_patch,land,convland]=land_gen2(matsize,patchsize,numpatch);
    elseif pturnover == 1 && patchsize == 1
        [patchID,n_patch,land]=land_gen(matsize,numpatch);
    elseif pturnover > 0 && pturnover < 1 && patchsize > 1
        [land,convland]=patch_turnover(pturnover,matsize,patchsize,convland,n_patch);
    elseif pturnover > 0 && pturnover < 1 && patchsize == 1
        [patchID,land]=patch_turnover1(matsize,patchID,n_patch,pturnover);
    end
    
    %% Germinate seeds
    fill_seed = rand(matsize);                              % random number matrix for determining seed germination
    SP(SP == 0 & (d1 > 0 | d2 > 0) & land==1 & ...          % species 2
        d2./(d1+d2) > fill_seed) = 2;
    SP(SP == 0 & (d1 > 0 | d2 > 0) & land==1 & ...          % species 1
        d2./(d1+d2) <= fill_seed) = 1;
    
    %% visualization
    if im_i==imtime
        figure(fig)
        subplot(1,2,1)
        imagesc(SP)
        colormap([1 1 1
            0 0.4470 0.7410
            0.8500 0.3250 0.0980]);
        cbh=colorbar('v');
        caxis([-0.5 2.5])
        set(cbh,'YTick',0:2)
        set(cbh,'YTickLabel',{'No Habitat','Species 1','Species 2'})
        im_i=0;
        
        subplot(1,2,2)
        plot(1:(i),cumsp1(1:(i)),1:(i),cumsp2(1:(i)),'LineWidth',2)
        axis([1,200,0,max(max(cumsp1),max(cumsp2))])
        xlabel('Time')
        ylabel('Population size')
        set(gca,'box','off'); %remove the box
        set(gca,'XTickLabel',[]) %remove the ticks
        set(gca,'YTickLabel',[]) %remove the ticks
        legend(strcat('\mu = ',num2str(mu(1))),strcat('\mu = ',num2str(mu(2))))
       
        drawnow
%         frame = getframe(fig);
%         im{idx} = frame2im(frame);       
%         idx=idx+1;
    end
    
    %% Determine whether simulations have converged
    if (sum(sum(SP==1)) < 0.05*sum(sum(land)) && sum(sum(SP==2))/sum(sum(land)) > 0.8) || ...
            (sum(sum(SP==2)) < 0.05*sum(sum(land)) && sum(sum(SP==1))/sum(sum(land)) > 0.8) % One species dominates
        
        convergence = 1;
        
    elseif i > 200
        if  ((sum(sum(SP==1)) < 0.05*sum(sum(land)) && ...                  % Populations are stable
                mean(cumsp2(i-200:i-100)) < mean(cumsp2(i-100:i))) || ...
                (sum(sum(SP==2)) < 0.05*sum(sum(land)) && ...
                mean(cumsp1(i-200:i-100)) < mean(cumsp1(i-100:i))))
            
            convergence = 1;
            
        end
    end
    
    im_i=im_i+1;
    i=i+1;
end

toc


% filename = 'testAnimated.gif'; % Specify the output file name
% for idx = 1:length(im)
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.2);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.2);
%     end
% end


