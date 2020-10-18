%%%
Plots behavior and simulated data of the ventriloquism effect. 


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
subj{1} = setdiff(1:24, [2 12]); % get rid of subjects who localized visual location
subj{2} = setdiff(1:24, [4 19 21]);


cols1 = [217,95,2 ; % YA
    27,158,119]./255;  % OA
    
ARG.VIS_scale = 8.5;
basecond = -4:4;

ms1 = 2;
ms = 5; % marker size


%% Get behavioral data
load(sprintf('%s/aging.mat',dat_dir),'data');
gvm = cell(length(grp), 1);
dva = unique(data{1}{1}(:, 4));
for g = 1:length(grp)
    for s = 1:length(data{g})
        for d = 1:length(dva)
            gvm{g}(s, d) = mean(data{g}{s}(data{g}{s}(:, 4)==dva(d), 6));
        end
    end
end

%% Load simulated data
ve = load(fullfile(dat_dir, 'model', 'VECI', 'VE_simSD.mat')); % brings xpsd{g, mdl}(b, 1, s): mean

%% Plot Figure
vis1 = figure;
vis1.Units = 'centimeters';
vis1.Position = [3 3 9 9];

fs = 12;
xlm1 = [-18 18];
xlm2 = [-35 35];


ofs = [4 4];
for g = 1:length(grp)
    subjs = subj{g};
    % Plot VE.
    lo = nan(1, length(basecond));
    hi = nan(1, length(basecond));
    tp = nan(1, length(basecond));
    for i = 1 : length(basecond)
        [lo(i), hi(i)] = confinth(gvm{g}(:, i), 'mean'); %, 0.05, 500);
        tp(i) = signrank(gvm{g}(:, i));        
    end
    p_corr = pval_adjust(tp, 'holm');
    p = patch([ARG.VIS_scale.*basecond fliplr(ARG.VIS_scale.*basecond)], [hi fliplr(lo)], cols1(g, :));
    set(p, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'Handlevisibility', 'off')
    hold on;
    plot(ARG.VIS_scale.*basecond.*ones(length(subjs), 1)+1.5-g, gvm{g}, 'color', cols1(g, :), 'LineStyle', 'none', 'Marker', 'o', ... 
    'MarkerSize', ms1, 'MarkerFaceColor', cols1(g, :), 'MarkerEdgeColor', 'none', 'HandleVisibility','off')
    hold on;
    xlim(xlm2)
    ylim([-34 34])
    a = ylim;

    plot(xlm2, [0 0], 'color',[0.6 0.6 0.6],'LineWidth',0.2, 'Handlevisibility', 'off')
    hold on;
    plot([0 0], a, 'color',[0.6 0.6 0.6], 'LineWidth',0.2, 'Handlevisibility', 'off')
    hold on;
    plot(ARG.VIS_scale.*basecond, nanmean(gvm{g}, 1), 'Color', cols1(g, :), 'LineWidth', 1.5, ...
        'LineStyle', '-', 'Marker', 's', ... 
        'MarkerSize', ms, 'MarkerFaceColor', cols1(g, :), 'MarkerEdgeColor', 'none', 'Handlevisibility', 'off')
    hold on; 
    simgvent = nanmean(ve.bxpsd{g,1}(:, 1, :), 3);
    plot(ARG.VIS_scale.*basecond, simgvent, 'Color', cols1(g, :), 'LineWidth', 1.5, ...
        'LineStyle', ':', 'Marker', 's', ... 
        'MarkerSize', ms, 'MarkerFaceColor', cols1(g, :), 'MarkerEdgeColor', 'none', 'Handlevisibility', 'off')

    set(gca, 'xtick', ARG.VIS_scale.*basecond(1:2:end), 'xticklabel', ARG.VIS_scale.*basecond(1:2:end), 'tickdir', 'out')

    box off
    hold all;

    if g==1
        plot(0, nan, 'color',cols1(1, :), 'LineWidth', 2)
        hold on; 
        plot(0, nan, 'color',cols1(2, :), 'LineWidth', 2)
        hold on;
        plot(0, nan, 'color',[0 0 0], 'LineStyle', '-','LineWidth', 2)
        hold on; 
        plot(0, nan, 'color',[0 0 0], 'LineStyle', ':','LineWidth', 2)
        legend('YA', 'OA','Real','Simulated','Location', 'NorthWest')
        legend boxoff

        xlabel('\DeltaVA (degree)', 'fontsize', fs)
        ylabel('ve (degree)')
    end
    offsetAxes(gca, ofs)
% Anne Urai (2020). offsetAxes(ax) (https://www.mathworks.com/matlabcentral/fileexchange/52351-offsetaxes-ax), 
% MATLAB Central File Exchange. Retrieved October 18, 2020.
end

filename = fullfile(fig_dir, 'Figure2');
print(vis1, filename, '-dtiff', '-r300')

