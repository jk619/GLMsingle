
path2fracridge = '/Users/jk7127/Documents/fracridge';
addpath(genpath(path2fracridge));
addpath(genpath('./../'))

if ~exist('./data','dir')
    mkdir('data')
    % Download TDM data for an example subject
    
    !curl -L --output ./data/information.mat https://osf.io/7hpr4/download
    !curl -L --output ./data/run1.mat https://osf.io/7jdzx/download
    !curl -L --output ./data/run2.mat https://osf.io/dhx3u/download
    !curl -L --output ./data/run3.mat https://osf.io/nrqz6/download
    !curl -L --output ./data/run4.mat https://osf.io/nwf28/download
    
end

%% load basic information and four runs of data
load ./data/information
numrun = dir('./data/run*');

data = cell(length(numrun),1);
for r = 1 : length(numrun)
    
    tmp = load(sprintf('./data/%s',numrun(r).name));
    data{r} = tmp.data;
    
end

%%
clc
disp(data)
fprintf('data consists of 4 runs, 870401 surface vertices, 6 cortical depths and 368 TRs for one example subject.\n\n')
disp(design(1:2)')
fprintf('Each run has a corresponding design matrix (TRs*conditions) with 6 different conditions \n\n')

figure(1);clf
subplot(2,2,1); imagesc(design{1}); colormap gray
xlabel('Conditions')
ylabel('TRs')
title('Design matrix for run1')
subplot(2,2,2); imagesc(design{2}); colormap gray
xlabel('Conditions')
ylabel('TRs')
title('Design matrix for run2')
subplot(2,2,3); imagesc(design{3}); colormap gray
xlabel('Conditions')
ylabel('TRs')
title('Design matrix for run3')
subplot(2,2,4); imagesc(design{4}); colormap gray
xlabel('Conditions')
ylabel('TRs')
title('Design matrix for run4')

%%
%We are going to split the data in segments of two runs. First analysis is
%going to be performed on run1 and run2 while the second analysis is going
%to be performed on run3 and run4. We do that to asses the reproduciblity
%of beta weights (run1+run2 vs. run3+run4) derived from GLM single. To speed up the estimaton we
%select data only from superficial layer and left V1.

load('./rois/lh_rois.mat');
V1_lh_ind = (lh_rois == 1 | lh_rois == 2); % dorsalV1 = 1 and ventralV1 = 2;

% create superficial dataset
superficial_data = cell(1,length(numrun)/2);
superficial_design = cell(1,length(numrun)/2);

for r = 1 : 2
    
    superficial_data{r} = squeeze(data{r}(V1_lh_ind,1,:));
    superficial_design{r} = design{r};
    
end
%%
% First, we run the GLM with canonical HRF only for run1+run2
opt.wantmemoryoutputs = [0 1 0 0]; % keep results in matlab memory 
opt.hrftoassume = 1;
opt.wantlibrary = 0;
opt.wantglmdenoise = 0;
opt.wantfracridge = 0;
outputidr = 'assume_HRF_test';
[assume_HRF_test] = GLMestimatesingletrial(superficial_design,superficial_data,stimdur,tr,outputidr,opt);

%%
% Now we enhance the model and include fitted HRF, GLM denoise and Ridge
% regression for run1+run2
opt.hrftoassume = 1; % get canonical hrf
opt.wantlibrary = 1; % create familly of HRFs output B
opt.wantglmdenoise = 1; % output C adds GLMDenoise
opt.wantfracridge = 1; % output D adds Ridge regression

opt.wantfileoutputs = [0 0 0 0]; % save results for output A, B, C, D
opt.wantmemoryoutputs = [0 1 1 1]; % keep results in matlab memory for A, B, C, D (note - We don't use output A as it's only ON/OFF GLM)

% run GLMsingle on first two runs
outputidr = 'fit_HRF_test';
[fit_HRF_test] = GLMestimatesingletrial(superficial_design,superficial_data,stimdur,tr,outputidr,opt);

%%
% for reproducibiliy we now run GLM on data from run3 and run4, first
% using canonical HRF only and later enhanciing with fitted HRF, GLMdenoise
% and Ridge regression

superficial_data = cell(1,length(numrun)/2);
superficial_design = cell(1,length(numrun)/2);
ct = 1;
for r = 3 : 4
    
    superficial_data{ct} = squeeze(data{r}(V1_lh_ind,1,:));
    superficial_design{ct} = design{r};
    ct = ct + 1;
end



%%
% We run the GLM with canonical HRF only for run3 + run4
opt.wantmemoryoutputs = [0 1 0 0]; % keep results in matlab memory for A, B, C, D
opt.hrftoassume = 1;
opt.wantlibrary = 0;
opt.wantglmdenoise = 0;
opt.wantfracridge = 0;
outputidr = 'assume_HRF_retest';
[assume_HRF_retest] = GLMestimatesingletrial(superficial_design,superficial_data,stimdur,tr,outputidr,opt);

%% We enhance the model to include fitted HRF, GLM denose and Ridge regression for run3 + run4
opt.hrftoassume = 1;
opt.wantlibrary = 1; % output B adds HRF fitting per voxel
opt.wantglmdenoise = 1; % output C adds GLMDenoie
opt.wantfracridge = 1; % output D adds Ridge regression

opt.wantfileoutputs = [0 0 0 0]; % save results for output A, B, C, D
opt.wantmemoryoutputs = [0 1 1 1]; % keep results in matlab memory for A, B, C, D
outputidr = 'fit_HRF_retest';
[fit_HRF_retest] = GLMestimatesingletrial(superficial_design,superficial_data,stimdur,tr,outputidr,opt);


%% Make voxel-vise correlation of beta weights for each approach (A,B,C,D)

approach = {'canonicalHRF';'fittedHRF';'GLMdenoise';'RR'};

figure(2);clf
for a = 1 : 4
    
    subplot(2,2,a)
    
    if a == 1
    plot(assume_HRF_test{2}.modelmd(:),assume_HRF_retest{2}.modelmd(:),'.')
    R = corr(assume_HRF_test{2}.modelmd(:),assume_HRF_retest{2}.modelmd(:));

    else
        plot(fit_HRF_test{a}.modelmd(:),fit_HRF_retest{a}.modelmd(:),'.')
        R = corr(fit_HRF_test{a}.modelmd(:),fit_HRF_retest{a}.modelmd(:));

    end
    xlabel('Beta test')
    ylabel('Beta retest')
    axis equal
    title(sprintf('R=%.2f, approach %s',R,approach{a}))
end
