% Load and preprocess
load('Scripts\Theta\Speed\theta_speed.mat');
t1.MidpointZT = (t1.StartZT + t1.EndZT) / 2;
t11 = t1(t1.state == "AWAKE" & t1.Speed > 1, :);

% Create time bins
edges = linspace(min(t11.MidpointZT), max(t11.MidpointZT), 4);
t11.TimeBin = discretize(t11.MidpointZT, edges, 'IncludedEdge','left');
t11.TimeBin = categorical(t11.TimeBin, 1:3, {'Bin1','Bin2','Bin3'});
t11.Condition = categorical(t11.Condition);
t11.Session = categorical(t11.Session);

% Fit mixed-effects model
formula = 'CentralFrequenyFooof ~ Speed * TimeBin * Condition + (1|Session)';
lme = fitlme(t11, formula);
disp(lme)
disp(anova(lme))

% Build prediction table
conds = categories(t11.Condition);
bins = categories(t11.TimeBin);
speed_ref = 0;
[groupCond, groupBin] = ndgrid(conds, bins);
pred_table = table;
pred_table.Condition = categorical(groupCond(:), conds);
pred_table.TimeBin = categorical(groupBin(:), bins);
pred_table.Speed = repmat(speed_ref, height(pred_table), 1);
pred_table.Session = repmat(t11.Session(1), height(pred_table), 1); % placeholder for random effect

% Predict with CI
[pred_table.PredictedFreq, pred_table.CI] = predict(lme, pred_table);
ci_half_width = diff(pred_table.CI, 1, 2) / 2;

% Plot predicted frequency
f=figure(2);clf; tiledlayout(1,2);f.Units="centimeters";
f.Position(3)=20;
f.Position(4)=6;
ax1=nexttile;hold on;ax2=nexttile;hold on;
colors = lines(numel(conds));
g = findgroups(pred_table.Condition);
axes(ax1); grid on;
for ic = 1:numel(conds)
    idx = g == ic;
    errorbar((1:3) + .2 * ic - .3, pred_table.PredictedFreq(idx), ci_half_width(idx), ...
        'Color', colors(ic,:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, 'DisplayName', conds{ic});
end
xticks(1:3); xticklabels(bins);
xlim([.5 3.5]);
grid on;
ylabel('Intercept (Theta Freq, Hz)');
xlabel('Time Bin');
legend show;
% title('Estimated Frequency Across Time Bins (Speed = 0 cm/s)');



% Extract fixed effects and confidence intervals
[beta, names, ci_fixed] = fixedEffects(lme, 'DFMethod', 'Residual');
names = cellstr(names.Name);  % Ensure 'names' is a cell array of character vectors

slopes = zeros(numel(conds), numel(bins));
slope_ci = zeros(numel(conds), numel(bins), 2);
get_idx = @(term) find(strcmp(names, term));
CovB = lme.CoefficientCovariance;

% Post-hoc bin comparisons
fprintf('\nPost-hoc Bin Comparisons Within Each Condition:\n');
for ic = 1:numel(conds)
    cond = conds{ic};
    fprintf('\nCondition: %s\n', cond);
    for i = 1:2
        for j = i+1:3
            bin_i = bins{i};
            bin_j = bins{j};

            new_data_i = table(categorical({cond}, conds), categorical({bin_i}, bins), speed_ref, ...
                t11.Session(1), 'VariableNames', {'Condition','TimeBin','Speed','Session'});
            new_data_j = new_data_i; new_data_j.TimeBin = categorical({bin_j}, bins);

            [pred_i, ci_i] = predict(lme, new_data_i);
            [pred_j, ci_j] = predict(lme, new_data_j);
            delta = pred_i - pred_j;

            se_diff = sqrt((diff(ci_i)/3.92)^2 + (diff(ci_j)/3.92)^2); % CI to SE approx
            tstat = delta / se_diff;
            pval = 2 * tcdf(-abs(tstat), lme.DFE);

            fprintf('  %s vs %s: Δ=%.3f Hz, p=%.4f\n', bin_i, bin_j, delta, pval);
        end
    end
end
%%
fprintf('\nPost-hoc Intercept Comparisons (SD vs NSD) Within Each Time Bin (FDR-corrected):\n');
pvals = zeros(numel(bins),1); deltas = zeros(numel(bins),1);
for ib = 1:numel(bins)
    bin = bins{ib};
    L = zeros(1, numel(beta));
    if bin == "Bin1"
        L(get_idx('Condition_SD')) = 1;
    elseif bin == "Bin2"
        L(get_idx('Condition_SD')) = 1;
        L(get_idx('Condition_SD:TimeBin_Bin2')) = 1;
    elseif bin == "Bin3"
        L(get_idx('Condition_SD')) = 1;
        L(get_idx('Condition_SD:TimeBin_Bin3')) = 1;
    end
    deltas(ib) = L * beta;
    se = sqrt(L * CovB * L');
    tval = deltas(ib) / se;
    pvals(ib) = 2 * tcdf(-abs(tval), lme.DFE);
end
[~, qvals] = fdr_bh(pvals);
for ib = 1:numel(bins)
    fprintf('  Bin %s: Δ=%.3f Hz (SD – NSD), p=%.4f, q=%.4f\n', bins{ib}, deltas(ib), pvals(ib), qvals(ib));
end
%%
fprintf('\nEstimated Slopes (Hz/cm/s) and 95%% CIs:\n');
for ic = 1:numel(conds)
    for ib = 1:numel(bins)
        cond = conds{ic};
        bin = bins{ib};

        % Construct contrast vector L
        L = zeros(1, numel(beta));
        idx_speed = get_idx('Speed');
        if isempty(idx_speed)
            error('Term "Speed" not found in fixed effects names. Found terms:\n%s', strjoin(names, '\n'));
        end

        L(get_idx('Speed')) = 1;
        if cond == "SD", L(get_idx('Condition_SD:Speed')) = 1; end
        if bin == "Bin2"
            L(get_idx('Speed:TimeBin_Bin2')) = 1;
            if cond == "SD"
                L(get_idx('Condition_SD:Speed:TimeBin_Bin2')) = 1;
            end
        elseif bin == "Bin3"
            L(get_idx('Speed:TimeBin_Bin3')) = 1;
            if cond == "SD"
                L(get_idx('Condition_SD:Speed:TimeBin_Bin3')) = 1;
            end
        end

        slope = L * beta;
        se_sl = sqrt(L * lme.CoefficientCovariance * L');  % <-- corrected
        tval = tinv(0.975, lme.DFE);  % 95% CI
        ci_L = [slope - tval * se_sl, slope + tval * se_sl];

        slopes(ic, ib) = slope;
        slope_ci(ic, ib, :) = ci_L;
        fprintf('  %s - %s: slope = %.3f Hz/cm/s (%.3f, %.3f)\n', ...
            cond, bin, slope, ci_L(1), ci_L(2));

    end
end

% Plot slopes with CIs
axes(ax2); grid on;
for ic = 1:numel(conds)
    x = (1:3) + .2 * ic - .3;
    y = slopes(ic, :);
    err = squeeze(diff(slope_ci(ic, :, :), 1, 3)) / 2;  % Transpose to 1×3

    if isvector(y) && isvector(err) && all(isfinite(y)) && all(isfinite(err))
        errorbar(x, y, err, 'Color', colors(ic,:), ...
            'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20, 'DisplayName', conds{ic});
    end
end
xticks(1:3); xticklabels(bins);grid on;xlim([.5 3.5]);

ylabel('Estimated Slope (Hz per cm/s)');
xlabel('Time Bin');
legend show;
% title('Speed–Theta Slope by Condition and Time Bin');
ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures/');
ff.save('speedLME')

%%
fprintf('\nPost-hoc Tests: Speed Slope ≠ 0 for Each Condition × TimeBin (FDR-corrected):\n');
pvals = zeros(numel(conds)*numel(bins),1);
row = 1;
for ic = 1:numel(conds)
    for ib = 1:numel(bins)
        cond = conds{ic}; bin = bins{ib};
        terms = {'Speed'};
        if cond == "SD", terms{end+1} = 'Condition_SD:Speed'; end
        if bin == "Bin2"
            terms{end+1} = 'Speed:TimeBin_Bin2';
            if cond == "SD", terms{end+1} = 'Condition_SD:Speed:TimeBin_Bin2'; end
        elseif bin == "Bin3"
            terms{end+1} = 'Speed:TimeBin_Bin3';
            if cond == "SD", terms{end+1} = 'Condition_SD:Speed:TimeBin_Bin3'; end
        end
        idx = cellfun(get_idx, terms);
        L = zeros(1, numel(beta)); L(idx) = 1;
        [F, pvals(row), ~] = coefTest(lme, L);
        row = row + 1;
    end
end
[~, qvals] = fdr_bh(pvals);
row = 1;
for ic = 1:numel(conds)
    for ib = 1:numel(bins)
        fprintf('  Condition: %s, Bin: %s → p=%.4f, q=%.4f\n', conds{ic}, bins{ib}, pvals(row), qvals(row));
        row = row + 1;
    end
end
%%

fprintf('\nCondition Differences in Speed Slope (SD vs NSD) by TimeBin (FDR-corrected):\n');
pvals = zeros(numel(bins),1); deltas = zeros(numel(bins),1);
for ib = 1:numel(bins)
    bin = bins{ib};
    terms = {'Condition_SD:Speed'};
    if bin == "Bin2", terms{end+1} = 'Condition_SD:Speed:TimeBin_Bin2'; end
    if bin == "Bin3", terms{end+1} = 'Condition_SD:Speed:TimeBin_Bin3'; end
    idx = cellfun(get_idx, terms);
    L = zeros(1, numel(beta)); L(idx) = 1;
    deltas(ib) = L * beta;
    se = sqrt(L * CovB * L');
    tval = deltas(ib) / se;
    pvals(ib) = 2 * tcdf(-abs(tval), lme.DFE);
end
[~, qvals] = fdr_bh(pvals);
for ib = 1:numel(bins)
    fprintf('  Bin: %s → Δ=%.3f Hz/cm/s, p=%.4f, q=%.4f\n', bins{ib}, deltas(ib), pvals(ib), qvals(ib));
end

%%
% Residual diagnostics
figure;
subplot(1,2,1);
plot(predict(lme), residuals(lme), 'k.');
xlabel('Fitted'); ylabel('Residuals');
title('Residuals vs Fitted'); grid on;

subplot(1,2,2);
qqplot(residuals(lme));
title('Normal Q-Q Plot');

%%

% Define range of speeds for prediction
speed_vals = linspace(1, 12, 50);  % from min to near max speed
nSpeeds = numel(speed_vals);

% Initialize prediction table
[condGrid, binGrid, speedGrid] = ndgrid(conds, bins, speed_vals);
pred_all = table;
pred_all.Condition = categorical(condGrid(:), conds);
pred_all.TimeBin = categorical(binGrid(:), bins);
pred_all.Speed = speedGrid(:);
pred_all.Session = repmat(t11.Session(1), height(pred_all), 1);  % any valid session

% Predict theta frequency
[pred_all.Freq, pred_all.CI] = predict(lme, pred_all);

% Plot
figure; hold on;
colors = lines(numel(conds));
linestyles = {'-','--','-.'};

for ic = 1:numel(conds)
    for ib = 1:numel(bins)
        idx = pred_all.Condition == conds(ic) & pred_all.TimeBin == bins(ib);
        sp = pred_all.Speed(idx);
        mu = pred_all.Freq(idx);
        ci = pred_all.CI(idx, :);

        % Plot with CI shading
        fill([sp; flipud(sp)], [ci(:,1); flipud(ci(:,2))], ...
            colors(ic,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        plot(sp, mu, linestyles{ib}, 'Color', colors(ic,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('%s - %s', conds{ic}, bins{ib}));
    end
end

xlabel('Speed (cm/s)');
ylabel('Predicted Theta Frequency (Hz)');
title('Theta Frequency vs Speed by Condition and TimeBin');
legend('Location','best'); grid on;


%%
fprintf('\nPost-hoc Pairwise Slope Comparisons (Hz/cm/s):\n');
names = cellstr(lme.CoefficientNames);  % ensure proper string array
get_idx = @(term) find(strcmp(names, term));
beta = fixedEffects(lme);
CovB = lme.CoefficientCovariance;

% Conditions and Bins
conds = categories(t11.Condition);
bins = categories(t11.TimeBin);

% --- 1. Bin vs Bin (within condition)
fprintf('\nWithin-condition Bin Slope Comparisons:\n');
for ic = 1:numel(conds)
    cond = conds{ic};
    fprintf('\nCondition: %s\n', cond);
    
    for i = 1:2
        for j = i+1:3
            bin_i = bins{i};
            bin_j = bins{j};

            % Contrast vectors for slope at bin i and j
            Li = zeros(1, numel(beta)); Lj = Li;
            Li(get_idx('Speed')) = 1;
            Lj(get_idx('Speed')) = 1;

            if cond == "SD"
                Li(get_idx('Condition_SD:Speed')) = 1;
                Lj(get_idx('Condition_SD:Speed')) = 1;
            end

            if bin_i == "Bin2"
                Li(get_idx('Speed:TimeBin_Bin2')) = 1;
                if cond == "SD", Li(get_idx('Condition_SD:Speed:TimeBin_Bin2')) = 1; end
            elseif bin_i == "Bin3"
                Li(get_idx('Speed:TimeBin_Bin3')) = 1;
                if cond == "SD", Li(get_idx('Condition_SD:Speed:TimeBin_Bin3')) = 1; end
            end

            if bin_j == "Bin2"
                Lj(get_idx('Speed:TimeBin_Bin2')) = 1;
                if cond == "SD", Lj(get_idx('Condition_SD:Speed:TimeBin_Bin2')) = 1; end
            elseif bin_j == "Bin3"
                Lj(get_idx('Speed:TimeBin_Bin3')) = 1;
                if cond == "SD", Lj(get_idx('Condition_SD:Speed:TimeBin_Bin3')) = 1; end
            end

            % Difference contrast
            Ldiff = Li - Lj;
            delta = Ldiff * beta;
            se = sqrt(Ldiff * CovB * Ldiff');
            tval = delta / se;
            pval = 2 * tcdf(-abs(tval), lme.DFE);
            fprintf('  %s vs %s: Δ=%.3f Hz/cm/s, p=%.4f\n', bin_i, bin_j, delta, pval);
        end
    end
end

% --- 2. Condition vs Condition (within bin)
fprintf('\nBetween-condition Slope Comparisons (SD vs NSD):\n');
for ib = 1:numel(bins)
    bin = bins{ib};

    % Contrast for SD minus NSD at given bin
    L = zeros(1, numel(beta));

    L(get_idx('Condition_SD:Speed')) = 1;  % SD - NSD slope difference

    if bin == "Bin2"
        L(get_idx('Condition_SD:Speed:TimeBin_Bin2')) = 1;
    elseif bin == "Bin3"
        L(get_idx('Condition_SD:Speed:TimeBin_Bin3')) = 1;
    end

    delta = L * beta;
    se = sqrt(L * CovB * L');
    tval = delta / se;
    pval = 2 * tcdf(-abs(tval), lme.DFE);
    fprintf('  Bin %s: Δ=%.3f Hz/cm/s (SD – NSD), p=%.4f\n', bin, delta, pval);
end


function [h, q] = fdr_bh(pvals, alpha)
    if nargin < 2, alpha = 0.05; end
    [p_sorted, sort_idx] = sort(pvals(:));
    m = length(p_sorted);
    thresh = (1:m)'/m * alpha;
    below = p_sorted <= thresh;
    max_id = find(below, 1, 'last');
    h = false(size(pvals)); q = zeros(size(pvals));
    if ~isempty(max_id)
        h(sort_idx(1:max_id)) = true;
    end
    q_sorted = p_sorted .* m ./ (1:m)';
    q_sorted = min(1, cummin(flip(q_sorted)));
    q(sort_idx) = flip(q_sorted);
end
