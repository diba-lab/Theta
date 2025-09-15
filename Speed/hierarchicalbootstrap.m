% Load data and calculate the midpoint of each bout
load('Scripts\Theta\Speed\theta_speed.mat');
%gets t11 from Fig4speed
t11=t1;annot='all';
t11=t1(t1.state=="AWAKE",:);annot='AWAKE';
t11=t1(t1.state=="AWAKE"&t1.Speed>1,:);annot='AWAKE, Speed > 1 cm/s';
% t11=t1(t1.Speed>1,:);annot='All, Speed > 1 cm/s';
hours(t11.StartZT)
t11.MidpointZT = (t11.StartZT + t11.EndZT) / 2;

% Time Binning
edges = linspace(min(t11.MidpointZT), max(t11.MidpointZT), 4);

% Initialize Hierarchical Bootstrap Analysis Parameters
n_iterations = 1000;
conditions = unique(t11.Condition);

% Initialize Tables for Storing Results
corr_results = array2table(zeros(3 * length(conditions), n_iterations), ...
    'VariableNames', "Corr_" + (1:n_iterations));
regression_results = cell(3 * length(conditions), 1);
combined_results = table(repelem(conditions, 3), repmat((1:3)', length(conditions), 1), ...
    corr_results, regression_results, 'VariableNames', {'Condition', 'TimeBin', 'Correlations', 'Regressions'});

% Hierarchical Bootstrap Analysis
for ibin = 1:length(edges)-1
    for icond = 1:length(conditions)
        % Filter data for this bin and condition
        condition_data = t11(t11.MidpointZT >= edges(ibin) & ...
            t11.MidpointZT < edges(ibin+1) & ...
            t11.Condition == conditions(icond), :);

        % Get list of sessions
        sessions = unique(condition_data.Session);

        bootstrap_correlations = zeros(n_iterations, 1);
        bootstrap_slopes = zeros(n_iterations, 1);
        bootstrap_intercepts = zeros(n_iterations, 1);

        for k = 1:n_iterations
            % Resample sessions with replacement
            sampled_sessions = sessions(randi(numel(sessions), numel(sessions), 1));
            sample = [];

            for s = sampled_sessions'
                bouts = condition_data(condition_data.Session == s, :);
                if isempty(bouts), continue; end
                % Resample bouts within session
                bout_idx = randi(height(bouts), height(bouts), 1);
                sample = [sample; bouts(bout_idx, :)];
            end

            % Proceed only if enough data
            if height(sample) >= 2
                [r, slope, intercept] = regressionAnalysis(sample.Speed, sample.CentralFrequenyFooof);
                bootstrap_correlations(k) = r;
                bootstrap_slopes(k) = slope;
                bootstrap_intercepts(k) = intercept;
            else
                bootstrap_correlations(k) = NaN;
                bootstrap_slopes(k) = NaN;
                bootstrap_intercepts(k) = NaN;
            end
        end

        % Store in table
        row_idx = (3 * (icond - 1) + ibin);
        combined_results.Correlations(row_idx, :) = array2table(bootstrap_correlations');
        combined_results.Regressions{row_idx} = table(bootstrap_slopes, bootstrap_intercepts, ...
            'VariableNames', {'Slope', 'Intercept'});
    end
end

%
% Confidence Interval Calculation
for ibin = 1:height(combined_results)
    % Confidence intervals for correlation coefficients
    ci_corr_values = prctile(table2array(combined_results.Correlations(ibin, :)), [2.5, 97.5]);
    combined_results.CI_Corr(ibin,:) = ci_corr_values;

    % Confidence intervals for slopes and intercepts
    reg_data = combined_results.Regressions{ibin};
    ci_slope_values = prctile(reg_data.Slope, [2.5, 97.5]);
    combined_results.CI_Slope(ibin,:) = ci_slope_values;

    ci_intercept_values = prctile(reg_data.Intercept, [2.5, 97.5]);
    combined_results.CI_Intercept(ibin,:) = ci_intercept_values;
end
% Statistical Test Between Bins 1-2, 2-3, and 1-3
p_values_corr = zeros(length(conditions), 3);
p_values_slope = zeros(length(conditions), 3);
p_values_intercept = zeros(length(conditions), 3);

for icond = 1:length(conditions)
    % Extract data for correlation coefficients
    bin1_corr = table2array(combined_results.Correlations((icond-1)*3 + 1, :));
    bin2_corr = table2array(combined_results.Correlations((icond-1)*3 + 2, :));
    bin3_corr = table2array(combined_results.Correlations((icond-1)*3 + 3, :));
    [~, p_corr_12] = ttest2(bin1_corr, bin2_corr); % Bin 1 vs Bin 2
    [~, p_corr_23] = ttest2(bin2_corr, bin3_corr); % Bin 2 vs Bin 3
    [~, p_corr_13] = ttest2(bin1_corr, bin3_corr); % Bin 1 vs Bin 3
    p_values_corr(icond, :) = [p_corr_12, p_corr_23, p_corr_13];

    % Extract data for regression slopes
    bin1_slope = combined_results.Regressions{(icond-1)*3 + 1}.Slope;
    bin2_slope = combined_results.Regressions{(icond-1)*3 + 2}.Slope;
    bin3_slope = combined_results.Regressions{(icond-1)*3 + 3}.Slope;
    [~, p_slope_12] = ttest2(bin1_slope, bin2_slope); % Bin 1 vs Bin 2
    [~, p_slope_23] = ttest2(bin2_slope, bin3_slope); % Bin 2 vs Bin 3
    [~, p_slope_13] = ttest2(bin1_slope, bin3_slope); % Bin 1 vs Bin 3
    p_values_slope(icond, :) = [p_slope_12, p_slope_23, p_slope_13];

    % Extract data for regression intercepts
    bin1_intercept = combined_results.Regressions{(icond-1)*3 + 1}.Intercept;
    bin2_intercept = combined_results.Regressions{(icond-1)*3 + 2}.Intercept;
    bin3_intercept = combined_results.Regressions{(icond-1)*3 + 3}.Intercept;
    [~, p_intercept_12] = ttest2(bin1_intercept, bin2_intercept); % Bin 1 vs Bin 2
    [~, p_intercept_23] = ttest2(bin2_intercept, bin3_intercept); % Bin 2 vs Bin 3
    [~, p_intercept_13] = ttest2(bin1_intercept, bin3_intercept); % Bin 1 vs Bin 3
    p_values_intercept(icond, :) = [p_intercept_12, p_intercept_23, p_intercept_13];
end

% Correct p-values for multiple comparisons using the Benjamini-Hochberg procedure
all_p_values = [p_values_corr(:); p_values_slope(:); p_values_intercept(:)];
% Perform Benjamini-Hochberg FDR correction
adj_p_values = benjaminiHochbergFDR(all_p_values);

% Helper function for Benjamini-Hochberg FDR correction
function adj_p_values = benjaminiHochbergFDR(p_values)
    % Sort p-values and get their original indices
    [sorted_p, sort_idx] = sort(p_values);
    n = length(p_values);
    % Calculate the Benjamini-Hochberg critical values
    bh_critical_values = (1:n)' / n * 0.05; % Assuming alpha = 0.05
    % Find the largest p-value that is less than the critical value
    below_threshold = sorted_p <= bh_critical_values;
    if any(below_threshold)
        max_idx = find(below_threshold, 1, 'last');
        threshold = sorted_p(max_idx);
    else
        threshold = 0;
    end
    % Adjust p-values
    adj_p_values = min(1, p_values * n / sum(p_values <= threshold));
end

% Reshape adjusted p-values back into their original matrices
adj_p_values_corr = reshape(adj_p_values(1:numel(p_values_corr)), size(p_values_corr));
adj_p_values_slope = reshape(adj_p_values(numel(p_values_corr)+1:numel(p_values_corr)+numel(p_values_slope)), size(p_values_slope));
adj_p_values_intercept = reshape(adj_p_values(numel(p_values_corr)+numel(p_values_slope)+1:end), size(p_values_intercept));

comparisons=array2table([p_values_corr, p_values_slope, p_values_intercept], ...
    'VariableNames', {'Corr_12', 'Corr_23', 'Corr_13', ...
                      'Slope_12', 'Slope_23', 'Slope_13', ...
                      'Intercept_12', 'Intercept_23', 'Intercept_13'}, ...
    'RowNames', cellstr(conditions));
%%
% Visualization
f=figure(2);clf; tiledlayout(1,3);f.Units="centimeters";
f.Position(3)=20;
f.Position(4)=6;
ax1=nexttile;hold on;ax2=nexttile;hold on;ax3=nexttile;hold on;


colors = colororder;
markers = {'o', '+', '*'};
for icond = 1:length(conditions)
    idx = (1:3) + (icond - 1) * 3;

    % Extract CI data for Correlation, Slope, and Intercept
    ci_corr = combined_results.CI_Corr(idx, :)';
    ci_slope = combined_results.CI_Slope(idx, :)';
    ci_intercept = combined_results.CI_Intercept(idx, :)';

    % Plot Correlation Coefficients
    axes(ax1); grid on;
    errorbar((1:3) + .2 * icond - .3, mean(ci_corr), (ci_corr(2, :) - ci_corr(1, :)) / 2, ...
        'Color', colors(icond, :), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20);

    % Add statistical significance markers for correlations
    sig_markers = {'Corr_12', 'Corr_23', 'Corr_13'};
    for i = 1:3
        if comparisons{char(conditions(icond)), sig_markers{i}} < 0.05
            if i == 3
                x_pos = 2 + .2 * icond - .3; % Above bin 2
            else
                x_pos = mean(i + [0, 1]) + .2 * icond - .3;
            end
            y_pos = mean(ci_corr(:)) + 0.05;
            text(x_pos, y_pos, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', colors(icond, :));
        end
    end

    % Plot Slopes
    axes(ax2); grid on;
    errorbar((1:3) + .2 * icond - .3, mean(ci_slope), (ci_slope(2, :) - ci_slope(1, :)) / 2, ...
        'Color', colors(icond, :), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20);

    % Add statistical significance markers for slopes
    sig_markers = {'Slope_12', 'Slope_23', 'Slope_13'};
    for i = 1:3
        if comparisons{char(conditions(icond)), sig_markers{i}} < 0.05
            if i == 3
                x_pos = 2 + .2 * icond - .3; % Above bin 2
            else
                x_pos = mean(i + [0, 1]) + .2 * icond - .3;
            end
            y_pos = mean(ci_slope(:)) + 0.05;
            text(x_pos, y_pos, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', colors(icond, :));
        end
    end

    % Plot Intercepts
    axes(ax3); grid on;
    errorbar((1:3) + .2 * icond - .3, mean(ci_intercept), (ci_intercept(2, :) - ci_intercept(1, :)) / 2, ...
        'Color', colors(icond, :), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20);

    % Add statistical significance markers for intercepts
    sig_markers = {'Intercept_12', 'Intercept_23', 'Intercept_13'};
    for i = 1:3
        if comparisons{char(conditions(icond)), sig_markers{i}} < 0.05
            if i == 3
                x_pos = 2 + .2 * icond - .3; % Above bin 2
            else
                x_pos = mean(i + [0, 1]) + .2 * icond - .3;
            end
            y_pos = mean(ci_intercept(:)) + 0.05;
            text(x_pos, y_pos, '*', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', colors(icond, :));
        end
    end
end

axes(ax1)
setAxes(combined_results, edges);
ylabel('R value');
axes(ax2)
setAxes(combined_results, edges);
ylabel('Slope');
axes(ax3)
setAxes(combined_results, edges);
ylabel('Intercept (Theta Freq, Hz)');

hold off;
ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures/');
annotation("textbox",[.1 .95 .8 .05],String=annot)
ff.save('speed')
%%

f=figure(1);clf;tiledlayout('horizontal');f.Units="centimeters";
f.Position(3)=20;
f.Position(4)=7;
ax(1)=nexttile;hold on;ax(2)=nexttile;hold on;ax(3)=nexttile;hold on;
t11.Duration = seconds(seconds(t11.EndZT-t11.StartZT));
t11.PowerFooofNorm = normalize(t11.PowerFooof,"range",[.1 1]);
colors=colororder;
for ibin = 1:length(edges)-1
    for icond = length(conditions):-1:1
        axes(ax(ibin));grid on;
        condition_data = t11(t11.MidpointZT >= edges(ibin) & t11.MidpointZT < edges(ibin+1) & ...
            t11.Condition == conditions(icond), :);
        s=scatter(condition_data.Speed, ...
            condition_data.CentralFrequenyFooof, ...
            seconds(condition_data.Duration),colors(icond,:),"filled");
        hold on;
        s.AlphaData=condition_data.PowerFooofNorm;
        s.AlphaDataMapping="none";
        s.MarkerFaceAlpha="flat";
        s.MarkerEdgeAlpha="flat";
        ax(ibin).XLim=[0 10];
        ax(ibin).XTick=0:2:20;
        ax(ibin).YLim=[5.6 8.4];
        [r, slope, intercept] = regressionAnalysis(condition_data.Speed, condition_data.CentralFrequenyFooof);
        pl=regressionPlot(condition_data.Speed, ...
            condition_data.CentralFrequenyFooof);
        pl.main.Color=colors(icond,:);
        pl.main.LineWidth=1.5;
        xlabel('Speed (cm/s)');ylabel('Frequency (Hz)');
        ax(ibin).DataAspectRatio=[1 range(ax(ibin).YLim)/range(ax(ibin).XLim) 1];
    end
    if ibin == 1
        % Add scatter points for legend
        legend_x = [8, 8, 8];
        legend_y = [6, 6.5, 7];
        legend_duration = [20, 50, 100];
        legend_colors = [0.5, 0.5, 0.5]; % Gray color for legend points

        axes(ax(ibin)); % Use the current axis
        s = scatter(legend_x, legend_y, legend_duration, legend_colors, "filled");
        s.MarkerFaceAlpha = 0.7; % Set transparency for legend points
    end
end
annotation("textbox",[.1 .90 .3 .1],String=annot)

ff.save('speed')
%% Take distros for speed<1

f=figure(3);clf;tiledlayout('horizontal');f.Units="centimeters";
t11=t1(t1.Speed<1,:);annot='All, Speed < 1 cm/s';
t11.MidpointZT = (t11.StartZT + t11.EndZT) / 2;

f.Position(3)=20;
f.Position(4)=8;
ax(1)=nexttile;hold on;ax(4)=nexttile;hold on;
ax(2)=nexttile;hold on;ax(5)=nexttile;hold on;
ax(3)=nexttile;hold on;ax(6)=nexttile;hold on;
t11.Duration = seconds(seconds(t11.EndZT-t11.StartZT));
t11.PowerFooofNorm = normalize(t11.PowerFooof,"range",[.1 1]);
colors=colororder;
for ibin = 1:length(edges)-1
    for icond = length(conditions):-1:1
        axes(ax(ibin));grid on;
        condition_data = t11(t11.MidpointZT >= edges(ibin) & t11.MidpointZT < edges(ibin+1) & ...
            t11.Condition == conditions(icond), :);
        s=scatter(condition_data.Speed, ...
            condition_data.CentralFrequenyFooof, ...
            seconds(condition_data.Duration),colors(icond,:),"filled");
        hold on;
        s.AlphaData=condition_data.PowerFooofNorm;
        s.AlphaDataMapping="none";
        s.MarkerFaceAlpha="flat";
        s.MarkerEdgeAlpha="flat";
        ax(ibin).XLim=[0 2];
        ax(ibin).XTick=0:2:20;
        ax(ibin).YLim=[5.6 8.4];
        [r, slope, intercept] = regressionAnalysis(condition_data.Speed, condition_data.CentralFrequenyFooof);
        pl=regressionPlot(condition_data.Speed, ...
            condition_data.CentralFrequenyFooof);
        pl.main.Color=colors(icond,:);
        pl.main.LineWidth=1.5;
        xlabel('Speed (cm/s)');ylabel('Frequency (Hz)');
        ax(ibin).DataAspectRatio=[1 range(ax(ibin).YLim)/range(ax(ibin).XLim)/5 1];

        axes(ax(ibin+3));grid on;
        histogram(condition_data.CentralFrequenyFooof,5:.25:10, ...
            FaceColor=colors(icond,:),EdgeColor=colors(icond,:));
        ax(ibin+3).XLim=[5.6 8.4];

        ax(ibin+3).View=[90 90];ax(ibin+3).XDir="reverse";
        ax(ibin+3).XTick=[];
        ylabel('# bouts');
    end
end
linkaxes(ax(4:6),'y')
annotation("textbox",[.1 .90 .3 .1],String=annot)

ff.save('speed distro s1')
%% plot single sessions
f=figure(3);clf;tiledlayout('horizontal');f.Units="centimeters";
t11=t1;annot='All';
t11=t1(t1.Speed>1,:);annot='All, Speed > 1 cm/s';
t11=t1(t1.state=="AWAKE"&t1.Speed>1,:);annot='AWAKE, Speed > 1 cm/s';
t11.MidpointZT = (t11.StartZT + t11.EndZT) / 2;

f.Position(3)=20;
f.Position(4)=8;
ax(1)=nexttile;hold on;
ax(2)=nexttile;hold on;
ax(3)=nexttile;hold on;
t11.Duration = seconds(seconds(t11.EndZT-t11.StartZT));
t11.PowerFooofNorm = normalize(t11.PowerFooof,"range",[.1 1]);
colors=colororder;
            sess=unique(t11.Session);
for ises=1:numel(sess)
    ses=sess(ises);
    t12=t11(t11.Session==ses,:);
    for ibin = 1:length(edges)-1
        for icond = length(conditions):-1:1
            axes(ax(ibin));grid on;
            condition_data = t12(t12.MidpointZT >= edges(ibin) & ...
                t12.MidpointZT < edges(ibin+1) & ...
                t12.Condition == conditions(icond), :);

            ses_data=condition_data(condition_data.Session==ses,:);
            if ~isempty(ses_data)
            s=scatter(ses_data.Speed, ...
                ses_data.CentralFrequenyFooof, ...
                seconds(ses_data.Duration),colors(icond,:),"filled");
            hold on;
            s.AlphaData=ses_data.PowerFooofNorm;
            s.AlphaDataMapping="none";
            s.MarkerFaceAlpha="flat";
            s.MarkerEdgeAlpha="flat";
            ax(ibin).XLim=[0 10];
            ax(ibin).XTick=0:2:20;
            ax(ibin).YLim=[5.6 8.4];
            [r, slope, intercept] = regressionAnalysis(ses_data.Speed, ...
                ses_data.CentralFrequenyFooof);
            pl=regressionPlot(ses_data.Speed, ...
                ses_data.CentralFrequenyFooof);
            pl.main.Color=colors(icond,:);
            pl.main.LineWidth=1.5;
            xlabel('Speed (cm/s)');ylabel('Frequency (Hz)');
            ax(ibin).DataAspectRatio=[1 range(ax(ibin).YLim)/...
                range(ax(ibin).XLim) 1];
            end
        end
    end
end
annotation("textbox",[.1 .90 .3 .1],String=annot)
ff.save('speed sessions')

% Helper function for linear regression analysis
function [r, slope, intercept] = regressionAnalysis(x, y)
r = corr(x, y);
coeffs = polyfit(x, y, 1);
slope = coeffs(1);
intercept = coeffs(2);
end
function [pl] = regressionPlot(x, y)
[x,I]=sort(x)
y=y(I);
r = corr(x, y);
[p,S] = polyfit(x, y, 1);
[y_fit,delta] = polyval(p,x,S);
% p(1)=plot(x,y,'bo');
pl.main=plot([x(1) x(end)],[y_fit(1) y_fit(end)],'r-');
% pl.ci=plot(x,y_fit+2*delta,'m--',x,y_fit-2*delta,'m--');
end

% Helper function to set Axes for plots
function setAxes(combined_results,edges)
ax = gca;
if isduration(edges), edges = hours(edges); end
ax.XTick = 1:3;
% ax.XTickLabel = arrayfun(@(x, y) sprintf('%.2f - %.2f', x, y), ...
%     edges(1:end-1), edges(2:end), 'UniformOutput', false);
ax.XLim = [.5, 3.5];
xlabel('Time Bin');
legend(unique(combined_results.Condition), 'Location', 'Best');
end
