clear;

mom_dir = fullfile('Z:', 'HP02');
dat_dir = fullfile(mom_dir, 'dat');
fig_dir = fullfile(mom_dir, 'fig');


addpath('C:\bootstrap') % requires the "Bootstrap Matlab Toolbox" 
%                                Ver 2.0
%
%               Abdelhak M. Zoubir and D. Robert Iskander

% --------------------------------------------------------------------
% PARAMETERS
%--------------------------------------------------------------------
grp = {'A','B'};
grpt = {'YA','OA'};
subj{1} = setdiff(1:24, [2 12]); % get rid of subjects who seemingly localized visual location
subj{2} = setdiff(1:24, [4 19 21]);

TruePosA = [-17,-8.5,0,8.5,17];
TruePosV = [-17,-8.5,0,8.5,17];
delta = [0 1 2 3 4];
vscale = 8.5;

models = {'ModAv','ModSel','Probmatch'};
modelst = {'MA','MS','PM'};

ARG = [];
ARG.ParinList = {'SA','SV1','SP','MP','PC'};
inc = [0.5 1 2];
i = 2;

cols1 = [217,95,2 ; % YA
    27,158,119]./255;  % OA

load(fullfile(dat_dir, 'nobs_A_AVtrial.mat')) % number of samples used to fit model per subject
load(fullfile(dat_dir, 'model','VECI','VE_BCI_result.mat')) % modeling result
%% Plot Figure

vis2 = figure;
vis2.Units = 'centimeters';
vis2.Position = [3 3 18 6];

fsa = 14;
annotation(vis2, 'textbox', [0.017, 0.99, 0.01, 0.01], 'String', '\bfA', 'FitBoxToText','on', 'LineStyle', 'none', 'fontsize', fsa)
annotation(vis2, 'textbox', [0.39, 0.99, 0.01, 0.01], 'String', '\bfB', 'FitBoxToText','on', 'LineStyle', 'none', 'fontsize', fsa)

fs = 9;

clear ax
% model frequency
ax(1) = subplot(1, 3, 1);
for g = 1:length(grp)
    bar((1:length(models))+0.1*(g-1), out{g}.Ef, 0.5, 'FaceColor', cols1(g, :), 'EdgeColor', 'none', 'ShowBaseLine', 'off')
    
    hold on;
end
ylim([0 1])
% grid on
xlim([0.5 length(models)+0.5])
plot(xlim, (1/length(models))*[1 1], '-g')
box off
legend(grpt)
legend boxoff
set(ax(1), 'xtick', 1:length(models), 'xticklabels', modelst, 'fontsize', fs,'fontweight','bold', 'tickdir', 'out')
ylabel('Model Frequency', 'fontsize', fs)
offsetAxes(ax(1), [4 4])
% Anne Urai (2020). offsetAxes(ax) (https://www.mathworks.com/matlabcentral/fileexchange/52351-offsetaxes-ax), 
% MATLAB Central File Exchange. Retrieved October 18, 2020.

% plot best model parameters
ylm = [-15 50];
ax(2) = subplot(1, 3, 2);
for g = 1:length(grp)    
    boxplot(bvp{g}(:, 1:length(ARG.ParinList)-1), 'positions',(1:length(ARG.ParinList)-1)+0.27*(g-1),'plotstyle','compact','colors', cols1(g, :), ...
        'medianstyle', 'target', 'outliersize', 0.1)
    set(gca, 'xticklabel', {' ',' ',' '}, 'activepositionproperty', 'position')
    hold on;

    for p = 1:length(ARG.ParinList)-1
        plot(0.27*(g-1)+p+0.1, bvp{g}(:, p), 'color', cols1(g, :), 'marker', 'o', 'markerfacecolor', cols1(g, :), ...
            'markeredgecolor', 'none', 'markersize', 3)
        hold on;
    end
end
hold on;
sx1 = [1.1 1.1; 1.1 1.37;1.37 1.37];
sy1 = [47 48; 48 48; 48 47];
sx2 = [2.1 2.1; 2.1 2.37; 2.37 2.37];
sy2 = [30 31; 31 31; 31 30];
plot(sx1, sy1, '-k'), hold on;
text(1.17, 49, '*', 'fontsize', 12), hold on;
plot(sx2, sy2, '-k'), hold on;
text(2.17, 32, '*', 'fontsize', 12)

set(ax(2), 'xtick', (1:length(ARG.ParinList)-1)+0.08, 'xticklabel', {'\sigma_A','\sigma_V','\sigma_P', '\mu_P'}, ...
    'fontsize', fs,'fontweight','bold', 'tickdir', 'out')
box off
xlim([0.6 length(ARG.ParinList)-1+0.6])
ylim(ylm)

ylabel('Position (degree)') 
set(ax(2), 'ytick',-10:10:50,'yticklabel', string(-10:10:50))


ax(3) = subplot(1, 3, 3);
for g = 1:length(grp)
    boxplot(bvp{g}(:, length(ARG.ParinList)), 'positions', 1+0.2*(g-1),'plotstyle','compact','colors', cols1(g, :), ...
        'medianstyle', 'target', 'outliersize', 0.1)
    set(ax(3), 'xticklabel', {' '}, 'activepositionproperty', 'position', 'YAxisLocation','right')
    hold on;
    p = length(ARG.ParinList);
    plot(0.2*(g-1)+1+0.08, bvp{g}(:, p), 'color', cols1(g, :), 'marker', 'o', 'markerfacecolor', cols1(g, :), ...
        'markeredgecolor', 'none', 'markersize', 3)
    hold on;
end
ylim([0 1])
xlim([0.6 1.6])

ylabel('Probability')
set(ax(3), 'xtick', 1+0.08, 'xticklabel', 'P_{COM}', 'fontsize', fs, 'tickdir', 'out')
box off
set(ax(3), 'ytick', 0:0.2:1, 'yticklabel', string(0:0.2:1))

ax(1).Position([1 2 3 4]) = [0.09 0.11 0.31 0.85];
ax(2).Position([1 2 3 4]) = [0.47 0.11 0.31 0.85]; % model params 1:4
ax(3).Position([1 2 3 4]) = [0.79 0.11 0.14 0.85]; % model params 5


filename = fullfile(fig_dir, 'Figure3');
print(gcf, filename, '-dtiff', '-r300')

