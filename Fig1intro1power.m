% Call the function to process and save Fooof data
            blocklist=categorical({'PRE','NSD','SD','TRACK','POST'});

processAndSaveFooofData(blocklist(:));

function processAndSaveFooofData(blocks)
    % Create an instance of SDFigures2 from the experiment.plot package
    sdf = experiment.plot.SDFigures2.instance;

    % Plot the Fooof data using the sdf instance
    % sdf.plotFooof;

    % Get theta peaks using the sdf instance
    tp = sdf.getThetaPeaks(blocks);

    % Get frequency time continuous data from theta peaks
    resc = tp.getFreqTimeContinuous;
    % Save the frequency time continuous data to a .mat file
    save('Scripts/Theta/resc.mat', 'resc');

    % Get frequency time bouts from theta peaks
    tblres = tp.getFreqTimebouts;
    % Save the frequency time bouts data to a .mat file
    save('Scripts/Theta/tblres.mat', 'tblres');
end

%%
% Call the function to process state ratios
processStateRatios;

function processStateRatios
    % Define the block category as "POST"
    blockcat = categorical("NSD");

    % Create an instance of SessionFactory from the experiment package
    sf = experiment.SessionFactory;

    % Define the selected sessions
    selected_ses = [1:2 4:12 14:15 18:23];

    % Get the session list for the selected sessions
    tses = experiment.SessionList(sf.getSessions(selected_ses), selected_ses);

    % Get the state ratios table for the selected sessions and block category
    stateratios = tses.getStateRatioTable(blockcat);

    % Save the state ratios table and session list to a .mat file
    save(sprintf('Scripts/Theta/stateratios-%s.mat', blockcat), 'stateratios', 'tses');
end
%%
% Load necessary data if not already in the workspace
dataFiles = {'tblres', 'resc', 'stateratios', 'tspeed', 'tses'};
filePaths = {'Scripts/Theta/tblres.mat', 'Scripts/Theta/resc.mat', ...
    'Scripts/Theta/stateratios-NSD.mat', 'Scripts/Theta/tspeed.mat', ...
    'Scripts/Theta/tses.mat'};

for i = 1:numel(dataFiles)
    if ~exist(dataFiles{i}, 'var')
        load(filePaths{i});
    end
end

% Clear unnecessary variables, keeping only the loaded data
clearvars -except tblres resc stateratios tses tspeed

% Initialize session factory and define states and conditions
sf = experiment.SessionFactory;
states = unique(tblres.state);
statesel = {'QWAKE', 'AWAKE'};
conds = unique(tblres.Condition);
condsel = {'NSD', 'SD'};

% Convert relevant columns to categorical
resc.Condition = categorical(resc.Condition);
resc.Block = categorical(resc.Block);
resc.StateTemp = categorical(resc.StateTemp);
resc.Session = categorical(resc.Session);

% Calculate awake ratio and add it to the stateratios table
awakeSum = sum([stateratios.("A-WAKE") stateratios.("Q-WAKE")], 2, "omitmissing");
sleepSum = sum([stateratios.("SWS") stateratios.("REM")], 2, "omitmissing");
awakeSum(isnan(awakeSum)) = 0; 
sleepSum(isnan(sleepSum)) = 0;
stateratios.awakeRatio = awakeSum ./ (awakeSum + sleepSum);

% Prepare figures and layouts
try close(f3); catch, end
f3 = figure(3); 
f3.Position = [2600 500 700 700];
tl3 = tiledlayout(3, 2);

try close(f1); catch, end
f1 = figure(1); 
f1.Position = [2600 500 700 550];
tl1 = tiledlayout(3, numel(condsel));

try close(f2); catch, end
f2 = figure(2); 
f2.Position = [2600 000 700 1000];
tl2 = tiledlayout(10, numel(condsel));

% Define colors and alpha values
Colors = orderedcolors("gem");
colorQA = Colors(4:5, :);
colorSD = Colors(1:2, :);
alphaQA = [0 1; .8 1.5];

% Define axis limits
ylim1 = [.2 1]; 
xlim1 = [0 5];
ylim2 = [.2 1];

% Normalize and filter tblres data
tblres.PowerFooofnorm = normalize(fillmissing(tblres.PowerFooof, ...
    "constant", mean(tblres.PowerFooof, "omitmissing")), "range", [.5 1]);
tblres.durnorm = normalize((tblres.end - tblres.start), "range", [10 70]);
tblres.durnorm = (tblres.end - tblres.start) / 2;
tblres(tblres.durnorm < 4, :) = [];
% Initialize axis index
axsi = 1;

for ic = 1:numel(condsel)
    tblcond = tblres(tblres.Condition == condsel{ic}, :);
    resccond = resc(resc.Condition == condsel{ic}, :);
    statecond = stateratios(stateratios.Condition == condsel{ic}, :);
    figure(f1);
    tf1(ic + numel(statesel)) = nexttile(tl1, ic + numel(statesel));
    tf1(ic) = nexttile(ic);

    for is = 1:numel(statesel)
        tblcondstate = tblcond(tblcond.state == statesel{is}, :);
        axes(tf1(ic));
        plotscatter(tblcondstate, colorQA(is, :), alphaQA(is, :));
        drawnow;

        sessions = unique(tblcondstate.Session);
        colors = colororder;

        for ises = 1:numel(sessions)
            tblcondstateses = tblcondstate(tblcondstate.Session == sessions{ises}, :);
            zt = arrayfun(@(ib) calculateZT(tblcondstateses(ib, :)), 1:height(tblcondstateses))';
            tblcondstateses = [tblcondstateses table(zt, 'VariableNames', {'ZT'})];

            resccondses = resccond(resccond.Session == sessions{ises}, :);
            sessionID = extractSessionID(sessions{ises});
            statecondses = statecond(statecond.Session == sessionID, :);
            figure(f2);
            tl = (ises - 1) * numel(condsel) + ic;
            tf2(tl) = nexttile(tl2, tl);
            scatterHandle = plotscatter(tblcondstateses, colorQA(is, :), alphaQA(is, :));
            ses = sf.getSessions(sessionID);

            if is == numel(statesel)
                tspeedsub = tspeed(tspeed.SessionID == sessionID, :);
                plotses(resccondses, [.3 .3 .3], 1:3, colorQA);
                ratiosubplots(tl) = plothyp(statecondses);

                if ises == 1
                    cb = colorbar("northoutside");
                    cb.Position(3) = cb.Position(3) / 3;
                    cb.Position(4) = cb.Position(4) * 2;
                    cb.Label.String = 'Awake Proportion';
                end

                axes(tf1(ic + numel(statesel)));
                plotses(resccondses, brighten(colors(ic, :), .3), 1, colorQA);
            end

            set(scatterHandle, 'ButtonDownFcn', @(src, event) ...
                updateAxes(src, event, tblcondstateses));
            if ises ~= numel(sessions)
                tf2(tl).XTickLabel = '';
            end
        end
    end

    [~,matp1]=plotsesmean(resccond, colors(ic, :), 1, colorQA);
    matp1.Time=hours(matp1.Time);
    st=statistics.fieldtrip.OneConditionComparison(matp1);
    vals1=mean(matp1.Data(:,1:6),2,"omitmissing");
            stat=st.getClusterBaseddepSamplesT(vals1);
        stat.report

    tlHype(ic) = plothypmean(statecond);
    hold off;
    figure(1);
    tf1(5) = nexttile(tl1, 5);
    plotsesmean(resccond, colors(ic, :), 1, colorQA); hold on;
    stat.plot;
    mean_vals = mean(vals1, 'omitnan');
    std_vals = std(vals1, 'omitnan') / sqrt(numel(vals1));thecolor=colors(ic,:);
    errorbar((.5), mean_vals, std_vals, 'Color', thecolor, 'LineWidth', 1.5); hold on;
    plot((.5), mean_vals, 'o', 'MarkerFaceColor', thecolor, 'MarkerSize', 6);
    tf1(6) = nexttile(tl1, 6);
    p{ic} = plotsesmean(resccond, colors(ic, :), 1:3, colorQA); hold on;
end

% Adjust axis labels and linking
adjustAxisLabels(tf1, tf2, tlHype, ratiosubplots, ylim1, xlim1, ylim2);

% Helper function to calculate ZT
function zt = calculateZT(row)
    try
        ti = row.EMG.getTimeInterval;
        zt = hours(row.startAbs - ti.getZeitgeberTime);
    catch
        zt = nan;
    end
end

% Helper function to extract session ID
function sessionID = extractSessionID(sessionStr)
    sescode = sessionStr; 
    sescode(1:3) = [];
    sessionID = str2double(sescode);
end

% Helper function to adjust axis labels and linking
function adjustAxisLabels(tf1, tf2, tlHype, ratiosubplots, ylim1, xlim1, ylim2)
    tf1(1).XTickLabel = '';
    tf1(2).XTickLabel = '';
    tf1(3).XTickLabel = '';
    tf1(4).XTickLabel = '';
    hold off;
    linkaxes(tf1);
    axes(tf1(1)); ylim(ylim1); xlim(xlim1);
    linkaxes([tf1 tlHype], 'x');
    linkaxes(tf2);
    axes(tf2(1)); ylim(ylim2); xlim(xlim1);
    linkaxes([tf2 ratiosubplots], 'x');
end
%%
try close(f4); catch, end; f4 = figure(4); f4.Position = [90   800   120   90];
plotscatterLegend(colorQA(1,:),alphaQA(1,:));
ax=gca;ax.Box="off";ax.YTick=[];ax.XTick=[];
try close(f5); catch, end; f5 = figure(5); f5.Position = [90   920   120   90];
plotscatterLegend(colorQA(2,:),alphaQA(2,:));
ax=gca;ax.Box="off";ax.YTick=[];ax.XTick=[];
%%
clearvars -except tblres resc stateratios tses
ff = logistics.FigureFactory.instance('C:\Users\ukaya\University of Michigan Dropbox\Utku Kaya\Kaya Sleep Project\Manuscript\Figures\Fig1-sup');
figure(1)
ff.save('Figure1S1power')
figure(2)
ff.save('Figure1power')
figure(4)
ff.save('Figure1L1power')
figure(5)
ff.save('Figure1L2power')

function [] = plotxcor(resc, stateratios)
    resc = sortrows(resc, {'Session', 'ZTCenter'});
    stateratios = sortrows(stateratios, {'Session', 'ZTCenter'});

    idx = resc.Condition == "NSD" & ~isnan(resc.CentralFrequency) & ...
        ~isnan(stateratios.awakeRatio);
    cf = resc(idx, :).CentralFrequency;
    ar = stateratios(idx, :).awakeRatio;
    [c, lags] = xcorr(cf, ar);
    nexttile(1);
    plot(lags * 3, c);
    xlabel('XCorr Lag (min)')
    xlim([-60 60]);
    nexttile(2)
    s = scatter(ar, cf, 'filled');
    s.MarkerFaceAlpha = .3;
    s.SizeData = 20;
    ylim([5 8.5]);
    xlabel('Awake Proportion');
    ylabel('Central Frequency (Hz)')

    a = pivot(resc, Columns = "Session", Rows = "ZTCenter", ...
        DataVariable = "CentralFrequency", IncludeMissingGroups = false);
    data = table2array(a(:, 2:end))';
    data(data == 0) = nan;
    resc1 = smoothdata(data, 2, 'gaussian', 10, 'omitmissing');

    a = pivot(stateratios, Columns = "Session", Rows = "ZTCenter", ...
        DataVariable = "awakeRatio", IncludeMissingGroups = false);
    data = table2array(a(:, 2:end))';
    data(data == 0) = nan;
    dist = unique(round(minutes(diff(a.ZTCenter)), 2));
    npoints = numel(a.ZTCenter);
    t = linspace(-dist * (npoints - 1), dist * (npoints - 1), 2 * npoints - 1);
    stateratios1 = smoothdata(data, 2, 'gaussian', 10, 'omitmissing');
    idxnan = any(isnan(stateratios1), 2) | any(isnan(resc1), 2);
    stateratios1(idxnan, :) = [];
    resc1(idxnan, :) = [];
    nexttile(3);
    y = 1:size(resc1);
    imagesc(hours(a.ZTCenter), y, resc1);
    xlabel('ZT Hours')
    ylabel('Sessions')
    xlim([0 5])
    nexttile(4)
    imagesc(hours(a.ZTCenter), y, stateratios1);
    xlabel('ZT Hours')
    ylabel('Sessions')
    nexttile(5)
    xlim([0 5])

    for i = 1:size(resc1, 1)
        [carr(i, :), lags] = xcorr(resc1(i, :), stateratios1(i, :));
    end

    plot(t, carr);
    xlabel('XCorr Lag (min)')
    xlim([-60 60]);

end

function [s] = plotscatter(tbl, color,lowcutoffforaplha)
    tblsorted = sortrows(tbl, "startAbs");
    zt = nan(height(tblsorted), 1);

    for ib = 1:height(tblsorted)
        ti = tblsorted(ib, :).EMG.getTimeInterval;
        try 
            zt(ib, 1) = hours(tblsorted(ib, :).startAbs - ti.getZeitgeberTime);
        catch ME
            zt(ib, 1) = nan;
        end
    end
    
    dur=tbl.end-tbl.start;
    s = scatter(zt, tbl.PowerFooof, dur*.5, "filled", ...
        AlphaData = (tbl.PowerFooof-lowcutoffforaplha(1))/diff(lowcutoffforaplha), ...
        MarkerFaceAlpha = "flat", ...
        AlphaDataMapping="none", ...
        MarkerFaceColor = color); hold on
end

function [s] = plotscatterLegend(color, lowcutoffforalpha)
% Define data
hold on;
x = [1, 2.5, 4, 1, 2.5, 4];
y = [10, 10, 10, 11, 11, 11];
duration = [50, 50, 50, 10, 50, 100];
m1=mean(lowcutoffforalpha);
pwrs1=linspace(lowcutoffforalpha(1)+.2,lowcutoffforalpha(2),3);
power = [pwrs1(1),pwrs1(2),pwrs1(3), m1, m1, m1];
labels = {num2str(pwrs1(1)), num2str(pwrs1(2)), num2str(pwrs1(3)), '20 s', '100 s', '200 s'}; % Legend labels

% Plot each point individually
for i = 1:length(x)
    s = scatter(x(i), y(i), duration(i), 'filled', ...
        'MarkerFaceColor', color, ...
        'MarkerFaceAlpha', (power(i)-lowcutoffforalpha(1))/diff(lowcutoffforalpha));
end

% Add legend
legend(labels, 'Location', 'best',NumColumns=2);

end
function [] = plotses(tbl, color, what, colorQA)
    list = {'RelativePower', 'RelativePowerQW', 'RelativePowerAW'};
    colors = [color; colorQA];
    mediansmooths = [10 30 30];
    linewidth = [2 .5 .5];

    for iline = what
        data = tbl.(list{iline});
        data = smoothdata(data, 1, 'movmedian', mediansmooths(iline), 'omitmissing');
        data = smoothdata(data, 1, 'gaussian', 10, 'omitmissing');
        p = plot(hours(tbl.ZTCenter), data); hold on
        p.LineWidth = linewidth(iline); p.Color = colors(iline, :);
    end

    drawnow;
end

function [] = plotspeed(tbl, color, colorQA)
    list = unique(tbl.States);
    colors = [color; colorQA];
    mediansmooths = [10 10 10];
    linewidth = [.5 .5 1.5];

    for iline = what
        tblsub = tbl(tbl.States == list(iline), :);
        data = tblsub.Speed;
        data = smoothdata(data, 1, 'movmedian', mediansmooths(iline), 'omitmissing');
        data = smoothdata(data, 1, 'gaussian', 10, 'omitmissing');
        p = plot(hours(tblsub.ZTCenter), data, ':'); hold on
        p.LineWidth = linewidth(iline); p.Color = colors(iline, :);
    end

    drawnow;
end

function [eb,matp1] = plotsesmean(tbl, color, what, colorQA)
    t = hours(unique(tbl.ZTCenter))';
    list = {'RelativePower', 'RelativePowerQW', 'RelativePowerAW'};
    colors = [color; colorQA];
    mediansmooths = [10 30 30];
    linewidth = [1 .5 .5];

    for iline = what
        a = pivot(tbl, Columns = "Session", Rows = "ZTCenter", ...
            DataVariable = list{iline}, IncludeMissingGroups = false);
        data1 = table2array(a(:, 2:end))';
        data1(data1 == 0) = nan;
        data1 = smoothdata(data1, 2, 'movmedian', mediansmooths(iline), 'omitmissing');
        data1 = smoothdata(data1, 2, 'gaussian', 10, 'omitmissing');
        err = std(data1, "omitmissing") / sqrt(size(data1, 1));
            matp1=data.basic.ChannelTime([],t,data1);

        eb(iline) = shadedErrorBar(t, mean(data1, "omitmissing"), err);
        eb(iline).mainLine.LineWidth = 2.5; eb(iline).mainLine.Color = colors(iline, :);
        eb(iline).patch.FaceColor = colors(iline, :); eb(iline).patch.EdgeColor = ...
            colors(iline, :);
        eb(iline).edge(1).Color = colors(iline, :);
        eb(iline).edge(2).Color = colors(iline, :);
    end

end

function ax1 = plothyp(tbl)
    ax = gca;
    ax1 = axes(Position = ax.Position);
    ax1.Position(2) = ax1.Position(2) + ax1.Position(4);
    ax1.Position(4) = ax1.Position(4) / 10;
    imagesc(hours(tbl.ZTCenter), 1, tbl.awakeRatio');
    cl = flipud(othercolor('RdBu9', 20));
    colormap(gca, cl)
    clim([0 1]);
    ax1.XAxis.Visible = "off";
    ax1.YAxis.Visible = "off";
    ax1.XTickLabel; drawnow;
    ax1.XLim=[0 5];
end

function ax1 = plothypmean(tbl)
    ax = gca;
    ax1 = axes(Position = ax.Position);
    ax1.Position(2) = ax1.Position(2) + ax1.Position(4);
    ax1.Position(4) = ax1.Position(4) / 10;
    a = pivot(tbl, Columns = "Session", Rows = "ZTCenter", ...
        DataVariable = "awakeRatio", IncludeMissingGroups = false);
    data = table2array(a(:, 2:end))';
    data(data == 0) = nan;
    data = smoothdata(data, 2, 'gaussian', 10, 'omitmissing');
    imagesc(hours(tbl.ZTCenter), 1, mean(data, "omitmissing"));
    cl = flipud(othercolor('RdBu9', 20));
    colormap(gca, cl)
    clim([0 1]);
    ax1.XAxis.Visible = "off";
    ax1.YAxis.Visible = "off";
    ax1.XTickLabel; drawnow;
        ax1.XLim=[0 5];

end
