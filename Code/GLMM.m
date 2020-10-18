%% 
clear 

%% Directories
mom_dir = fullfile('Z:', 'HP02');
log_dir = fullfile(mom_dir, 'log');
dat_dir = fullfile(mom_dir, 'dat');
fig_dir = fullfile(mom_dir, 'fig');

addpath('Y:\Matlab\ckmatlab\ckstatistics')
addpath('Z:\HP02\VBA-toolbox')

%% Single subject regression
grp = {'A','B'};
grpt = {'YA','OA'};
subj{1} = setdiff(1:24, [2 12]);
subj{2} = setdiff(1:24, [4 19 21]);

load(sprintf('%s/aging.mat',dat_dir),'data');

% data
% c1: AV trial V stimulus location
% c2: AV trial A stimulus location
% c3: A trial A stimulus location
% c4: dVA: AV V - AV A (a.k.a, audio-visual discrepancy)
% c5: AV trial response
% c6: VE
% c7: A trial response
% c8: VAE

%% compare effects of previous dVA (combined data for both groups)
clear gdata tdata ttdata
gdata = [];
for g=1:length(data)
    ttdata = [];
    for s=1:length(data{g})
        tdata = data{g}{s}; 
        tdata = cat(2, sign(tdata(:, 4)).*sqrt(abs(tdata(:, 4))), tdata(:, 4), sign(tdata(:, 5)).*sqrt(abs(tdata(:, 5))), ... 
            tdata(:, 5), tdata(:, [6 8]), g*ones(size(tdata, 1), 1), (g-1)*length(data{1})+s*ones(size(tdata, 1), 1)); 
        ttdata = cat(1, ttdata, tdata);
    end
    gdata = cat(1, gdata, ttdata);
end
clear Table
Table = array2table(gdata, 'VariableNames', {'qDp', 'Dp', 'qRp','Rp','VE', 'VAE', 'grp', 'subj'});
Table.grp = categorical(Table.grp);

clear GModel3
% VE eff
GModel3{1,1} = fitglme(Table, ...
    'VE ~ 1 + qDp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
GModel3{1,2} = fitglme(Table, ...
    'VE ~ 1 + Dp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 
GModel3{1,3} = fitglme(Table, ...
    'VE ~ 1 + qDp + Dp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 

% VAE eff
GModel3{2,1} = fitglme(Table, ...
    'VAE ~ 1 + qDp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
GModel3{2,2} = fitglme(Table, ...
    'VAE ~ 1 + Dp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 
GModel3{2,3} = fitglme(Table, ...
    'VAE ~ 1 + qDp + Dp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 
GModel3{2,4} = fitglme(Table, ...
    'VAE ~ 1 + Dp + Rp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 
GModel3{2,5} = fitglme(Table, ...
    'VAE ~ 1 + Rp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 

BIC = [];
for eff = 1:2
    for nl = 1:3
        BIC(eff, nl) = GModel3{eff, nl}.ModelCriterion.BIC;
    end
end
round(BIC-min(BIC, [], 2))

%    108   112     0
%     10     0     9

% now for VAE, compare with model with Rp
BIC = [];
for eff = 2
    for nl = 1:4
        BIC(eff, nl) = GModel3{eff, nl}.ModelCriterion.BIC;
    end
end
round(BIC(2, [2 4])-min(BIC(2, [2 4]), [], 2))
% 74     0 -> VAE ~ 1 + Dp + Rp + (1|subj) is the best.
% BF
log10(exp((GModel3{2,2}.ModelCriterion.BIC-GModel3{2,4}.ModelCriterion.BIC)/2))
% 16.0101

%% with fitglme, single trial regression => see if VE and VAE differs between groups
clear gdata tdata ttdata
gdata = [];
for g=1:length(data)
    ttdata = [];
    for s=1:length(data{g})
        tdata = data{g}{s}; 
        tdata = cat(2, sign(tdata(:, 4)).*sqrt(abs(tdata(:, 4))), tdata(:, 4), sign(tdata(:, 5)).*sqrt(abs(tdata(:, 5))), ... 
            tdata(:, 5), tdata(:, [6 8]), g*ones(size(tdata, 1), 1), (g-1)*length(data{1})+s*ones(size(tdata, 1), 1)); 
        ttdata = cat(1, ttdata, tdata);
    end
    gdata = cat(1, gdata, ttdata);
end
clear Table
Table = array2table(gdata, 'VariableNames', {'qDp', 'Dp', 'qRp','Rp','VE', 'VAE', 'grp', 'subj'});
Table.grp = categorical(Table.grp);

clear GModel1
% VE eff
GModel1{1,1} = fitglme(Table, ...
    'VE ~ 1 + qDp + Dp + qDp:grp + Dp:grp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
% GModel1{1,2} = fitglme(Table, ...
%     'VE ~ 1 + qDp*grp + Dp*grp + Rp*grp + (1|subj)', ...
%     'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
%     'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); % -> this model is regressing against itself (VE). 

% VAE eff
GModel1{2,1} = fitglme(Table, ...
    'VAE ~ 1 + Dp + Dp:grp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 
GModel1{2,2} = fitglme(Table, ...
    'VAE ~ 1 + Dp + Rp + Dp:grp + Rp:grp + (1|subj)', ...
    'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
    'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 


BIC = [];
for eff = 2
    for nl = 1:2
        BIC(eff, nl) = GModel1{eff, nl}.ModelCriterion.BIC;
    end
end
BIC-min(BIC, [], 2)
% vae 93.8316         0


%% get slopes for groups sepearately 
clear GModel2
for g=1:length(data)
    ttdata = [];
    for s=1:length(data{g})
        tdata = data{g}{s}; % we need c4, c6, c8
        tdata = cat(2, sign(tdata(:, 4)).*sqrt(abs(tdata(:, 4))), tdata(:, 4), sign(tdata(:, 5)).*sqrt(abs(tdata(:, 5))), ... 
            tdata(:, 5), tdata(:, [6 8]), s*ones(size(tdata, 1), 1));
        ttdata = cat(1, ttdata, tdata);
    end
    clear Table
    Table = array2table(ttdata, 'VariableNames', {'qDp', 'Dp', 'qRp', 'Rp', 'VE', 'VAE', 'subj'});

    % VE eff
    GModel2{g}{1,1} = fitglme(Table, ...
        'VE ~ 1 + qDp + (1|subj)', ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
        'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
    GModel2{g}{1,2} = fitglme(Table, ...
        'VE ~ 1 + Dp + (1|subj)', ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
        'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
    GModel2{g}{1,3} = fitglme(Table, ...
        'VE ~ 1 + qDp + Dp + (1|subj)', ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
        'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);

    % VAE eff
    GModel2{g}{2,1} = fitglme(Table, ...
        'VAE ~ 1 + Dp + (1|subj)', ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
        'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
    GModel2{g}{2,2} = fitglme(Table, ...
        'VAE ~ 1 + Rp + (1|subj)', ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
        'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20);
    GModel2{g}{2,3} = fitglme(Table, ...
        'VAE ~ 1 + Dp + Rp + (1|subj)', ...
        'Distribution', 'Normal', 'Link', 'Identity', 'FitMethod', 'Laplace',...
        'DummyVarCoding', 'reference', 'PLIterations', 500, 'InitPLIterations', 20); 
end

BIC = {};
for g = 1:2
    for eff = 1:2
        for nl = 1:3
            BIC{g}(eff, nl) = GModel2{g}{eff, nl}.ModelCriterion.BIC;
        end
    end
end
round(BIC{1}-min(BIC{1},[], 2))
%     1.9290  105.6365         0 % ve
%    19.7694  538.5363         0 % vae
round(BIC{2}-min(BIC{2},[], 2))
%   129.8679   19.9378         0 % ve
%    65.3693   96.0384         0 % vae


%% Load GLMM data and calculate BF
% load(fullfile(dat_dir, 'GModel.mat'))
%%%%%%%%%%%%% VE %%%%%%%%%%%%%%%%%
BF10_FM = [];
BF10_FM(1) = log10(exp((GModel3{1,2}.ModelCriterion.BIC-GModel1{1,1}.ModelCriterion.BIC)/2)); % BF for qDp
BF10_FM(2) = log10(exp((GModel3{1,1}.ModelCriterion.BIC-GModel1{1,1}.ModelCriterion.BIC)/2)); % BF for Dp
BF10_FM(3) = log10(exp((GModel3{1,3}.ModelCriterion.BIC-GModel1{1,1}.ModelCriterion.BIC)/2)); % BF for G
% 114.5649  113.7283   90.2502

% for each group, efficacy of qDp and Dp
BF10_GL = [];
for g = 1:2
    BF10_GL(g, 1) = log10(exp((GModel2{g}{1,1}.ModelCriterion.BIC-GModel2{g}{1,3}.ModelCriterion.BIC)/2)); % BF for Dp
    BF10_GL(g, 2) = log10(exp((GModel2{g}{1,2}.ModelCriterion.BIC-GModel2{g}{1,3}.ModelCriterion.BIC)/2)); % BF for qDp
    BF10_GL(g, 3) = log10(exp((GModel2{g}{1,2}.ModelCriterion.BIC-GModel2{g}{1,1}.ModelCriterion.BIC)/2)); % Dp vs. qDp
end
%     0.4189   22.9387   22.5198
%    28.2005    4.3294  -23.8710

%%%%%%%%%%%%% VAE %%%%%%%%%%%%%%%%%
clear BF10_FM
% full model: no Rp vs. Rp -> VAE ~ 1 + Dp + (1 | subj) vs. VAE ~ 1 + Dp + Rp + (1 | subj) 
% full model: no group vs. group -> VAE ~ 1 + Dp + Rp + (1 | subj) vs. VAE ~ 1 + Dp + Rp + Dp:grp + Rp:grp + (1 | subj) 
BF10_FM(1) = log10(exp((GModel3{2,5}.ModelCriterion.BIC-GModel1{2,2}.ModelCriterion.BIC)/2)); % BF for Dp
BF10_FM(2) = log10(exp((GModel3{2,2}.ModelCriterion.BIC-GModel1{2,2}.ModelCriterion.BIC)/2)); % BF for Rp
BF10_FM(3) = log10(exp((GModel3{2,4}.ModelCriterion.BIC-GModel1{2,2}.ModelCriterion.BIC)/2)); % BF for G
% BF10_FM =
% 
%    128.2881   26.7472   10.7371

% for each group, efficacy of Dp and Rp
BF10_GL = [];
for g = 1:2
    BF10_GL(g, 1) = log10(exp((GModel2{g}{2,2}.ModelCriterion.BIC-GModel2{g}{2,3}.ModelCriterion.BIC)/2)); % BF for Dp
    BF10_GL(g, 2) = log10(exp((GModel2{g}{2,1}.ModelCriterion.BIC-GModel2{g}{2,3}.ModelCriterion.BIC)/2)); % BF for Rp
    BF10_GL(g, 3) = log10(exp((GModel2{g}{2,2}.ModelCriterion.BIC-GModel2{g}{2,1}.ModelCriterion.BIC)/2)); % Dp vs. Rp
end
%   116.9417    4.2929  112.6488
%    20.8545   14.1948    6.6597