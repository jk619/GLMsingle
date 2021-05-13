
%% addpath to GLMsingle and fracridge
addpath(genpath('/Users/jankurzawski/Documents/GLMsingle'))
addpath(genpath('/Users/jankurzawski/Documents/fracridge'))

clear

if ~exist('./data','dir')
    mkdir('data')
    % Download data for an example subject
    
    !curl -L --output ./data/nsdflocexampledataset.mat https://osf.io/8rmjk/
    !curl -L --output ./data/haxby.mat https://osf.io/f6vpq/download
    !curl -L --output ./data/nsdcoreexampledataset.mat https://osf.io/jqpwz/download
end

datasets = {'nsdcore';'nsdfloc';'haxby'};
dataset = datasets{1};


switch dataset
    
    case 'nsdfloc'
        load('./data/nsdflocexampledataset.mat')
    case 'nsdcore'
        load('./data/nsdcoreexampledataset.mat')
    case 'haxby'
        load('./data/haxby.mat')
end

%%
clc
whos
% data -> Consists of several runs of 4D volume files (x y z t)  where 
% (t)ime is the 4th dimention

% design -> Each run has a corresponding design matrix where each colum
% describes single condition (conditions are repeated across runs). Each
% design matrix is binary with 1 specfing the time (TR) when stimulus is
% presented on the screen.


fprintf('There are %d runs in total.\n',length(design));
fprintf('The dimensions of the data for the first run are %s.\n',mat2str(size(data{1})));
fprintf('The stimulus duration is %.6f seconds.\n',stimdur);
fprintf('The sampling rate (TR) is %.6f seconds.\n',tr);

%%
figure(1);clf
for d = 1 : 4
    subplot(2,2,d); imagesc(newdesign{d}); colormap gray
    xlabel('Conditions')
    ylabel('TRs')
    title(sprintf('Design matrix for run%i',d))
    if exist('stimuli','var')
        xticks([1:length(stimuli)])
        xticklabels(stimuli)
        %
        
    end
    set(gca,'YDir','normal')
    
end
%%
figure(2);clf
imagesc(makeimagestack(data{1}(:,:,:,1)));
colormap(gray);
axis equal tight;
colorbar;
title('fMRI data (first volume)');
%% Call GLMestimatesingletrial using default parameters
% Outputs will be stored in 'results' and Figures will be written to
% 'example1figures' in the current directory.

% default options (all parameters below can be assigned to a variable i.e
% opt by creating fields (opt.wantlibrary). Options are the 6th input to
% GLMestimatesingletrial.

% wantlibrary = 1 -> Fit hRF to each voxel
% wantglmdenoise = 1 -> Use GLMdenoise
% wantfracridge = 1  -> Use ridge regression to improve beta estimates
% chunknum = 5000 -> is the number of voxels that we will process at the
% same time. For setups with lower memory deacrease this number.

% wantmemoryoutputs is a logical vector [A B C D] indicating which of the
%     four model types to return in the output <results>. The user must be careful with this,
%     as large datasets can require a lot of RAM. If you do not request the various model types,
%     they will be cleared from memory (but still potentially saved to disk).
%     Default: [0 0 0 1] which means return only the final type-D model.

% wantfileoutputs is a logical vector [A B C D] indicating which of the
%     four model types to save to disk (assuming that they are computed).
%     A = 0/1 for saving the results of the ONOFF model
%     B = 0/1 for saving the results of the FITHRF model
%     C = 0/1 for saving the results of the FITHRF_GLMDENOISE model
%     D = 0/1 for saving the results of the FITHRF_GLMDENOISE_RR model
%     Default: [1 1 1 1] which means save all computed results to disk.


[results] = GLMestimatesingletrial(design,data,stimdur,tr,'example1figures');
%% Important outputs
% <R2> is model accuracy expressed in terms of R^2 (percentage).
% <modelmd> is the full set of single-trial beta weights (X x Y x Z x TRIALS)
% <HRFindex> is the 1-index of the best HRF, can be retried with
% getcanonicalhrflibrary(stimdur,tr)
% <FRACvalue> is the fractional ridge regression regularization level chosen for each voxel
slice = 1;
val2plot = {'meanvol';'R2';'HRFindex';'FRACvalue'};
cmaps = {gray;hot;parula;copper};
figure(3);clf

for v = 1 : length(val2plot)
    
    f=subplot(2,2,v);
    imagesc(results{4}.(val2plot{v})(:,:,slice)); axis off image;
    colormap(f,cmaps{v}) % Error message is related to this line
    colorbar
    title(val2plot{v})
    set(gca,'FontSize',20)
    
end

set(gcf,'Position',[ 1224         840         758         408])

