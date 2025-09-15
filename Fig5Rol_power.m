

%% plot running Windows
blockstr=categorical("PRE");xlim1=[-4 0];
% blockstr=categorical("RUN");xlim1=[0 1];
% blockstr=categorical("NSD");xlim1=[0 5];
% blockstr=categorical("RS");xlim1=[0 5];

injlist=categorical({'CTRL','ROL'});
    states_sub=categorical({'AWAKE','QWAKE'});

conditionlist=categorical({'SD'});
s=load(sprintf('Scripts/Theta/matfiles-rol/runningWindowTheta-%s-%s.mat', ...
    join(string(states_sub)),join(string(blockstr))));
runningWindowTheta=s.runningWindowTheta;
s=load(sprintf('Scripts/Theta/matfiles-rol/stateratios-%s.mat', ...
    join(string(blockstr))));
stateratios=s.stateratios;
tses=s.tses;
injections = arrayfun(@(s) s.SessionInfo.Injection, tses, 'UniformOutput', false);
conditions = arrayfun(@(s) s.SessionInfo.Condition, tses, 'UniformOutput', false);
ids = (1:numel(tses))';  % generate indices
st = table(categorical(injections(:)), categorical(conditions(:)), ids, 'VariableNames', {'INJECTION', 'SLEEP', 'ID'});
ctrlidx=find(st.INJECTION==injlist(1)&st.SLEEP==conditionlist(1));
rolidx=find(st.INJECTION==injlist(2)&st.SLEEP==conditionlist(1));
% runningWindowTheta=runningWindowTheta(runningWindowTheta.Condition=="NSD",:);
f=figure(1);clf;f.Units="centimeters"; f.Position(3)=20;f.Position(4)=17;
tl=tiledlayout('vertical');
panel1=nexttile;grid on;hold on;
panel2=nexttile; grid on;hold on;

% axes(panel1)
% xlim(xlim1)
% ylim([5.5 8.5]);
% panel1.DataAspectRatio=[2 1 1];
% plothypmean(stateratios(ismember(stateratios.SessionIdx,find(ctrlidx)),:));
% xlim(xlim1);
% axes(panel2);
% xlim(xlim1);
% ylim([5.5 8.5]);
% panel2.DataAspectRatio=[2 1 1];
% plothypmean(stateratios(ismember(stateratios.SessionIdx,find(rolidx)),:));
% xlim(xlim1);
% ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures/f5');
% ff.save(sprintf('AwakeProportionBars-%s',blockstr))

    conditions=sort(unique(runningWindowTheta.Condition),'descend');
sessions=unique(runningWindowTheta.Session);
colors=colororder;
if blockstr=="RS"
    x=hours(round((0:1/20:5.5)*20)/20)';
else
    x=(min(runningWindowTheta.ZTCenter):hours(1/20):max( ...
        runningWindowTheta.ZTCenter))';
end
cfs=nan(numel(sessions),numel(x));
for icond=1:numel(conditionlist)
    cond=conditions(icond);
    cond_idx=st.SLEEP==cond;
    st_sub=st(cond_idx,:);
    injs=unique(st_sub.INJECTION);
    for iinj=1:numel(injlist)
        injselected=injlist(iinj);
        st1_sub=st_sub(st_sub.INJECTION==injselected,:);
        sub_ses=st1_sub.ID;
        axes(panel1);
        for ises=1:numel(sub_ses)
            sesno=sub_ses(ises);
            if blockstr=="NSD"||blockstr=="RUN"||blockstr=="RS"
                ses=tses(sesno);
                blockTimeWin=hours(ses.getBlockZT(blockstr));
                ztcorrection=hours(round(blockTimeWin(1)*20)/20);
            else
                ztcorrection=hours(0);
            end
            tbl_ses=runningWindowTheta(runningWindowTheta.Session==sesno,:);
            cf=tbl_ses.power;
            ztCenterCorrected=hours(round(hours(( ...
                tbl_ses.ZTCenter-ztcorrection)*20))/20);
            stateratios(ismember(stateratios.SessionIdx,sesno),:).ZTCenter=...
                hours(round(hours(( ...
                stateratios(ismember(stateratios.SessionIdx,sesno),:).ZTCenter -ztcorrection)*20))/20);
            hold on;
            idx=ismember(hours(x),hours(ztCenterCorrected));
            cfs(ismember(sessions,sesno),idx)=cf';
        end
        colorCond=colors(injlist==injselected,:);
        sub_cfs=cfs(ismember(sessions,sub_ses),:);
        matp1=data.basic.ChannelTime([],x,sub_cfs);
        matp2=matp1.getMedianFiltered(minutes(30));
        ctdatas{iinj}=matp2.getGaussianFiltered(minutes(30));
        axes(panel1);
        p=ctdatas{iinj}.plotChannels(colors(iinj,:));
        axes(panel2);
        p=ctdatas{iinj}.plotErrorBar(colors(iinj,:));
        legendObjects(iinj)=p.mainLine;
        xlim(xlim1);
        panel1.DataAspectRatio=[5 1 1];
        panel2.DataAspectRatio=[5 1 1];
        panel1.XTick=-5:10;
        panel2.XTick=-5:10;
        panel2.YLim=[0.0 1.0];
        linkaxes([panel1 panel2],'xy')
        if blockstr=="RUN"
            xlim1=hours([0 max(x)]);
        end
        xlim(xlim1);
    end
end
xlabel('ZT (h)');ylabel('Relative Power (mV²/Hz)')
if blockstr=='NSD'
    image = imread('Scripts/Theta/injection.jpg');
    injtimes=[1 2.5];
    ax=gca;
    vline(injtimes)
    % imshow(ax,image,'XData', [injtimes(1), injtimes(1)+.5 ], 'YData', [min(ax.YLim) min(ax.YLim)+0.5])
    % axis on
end
axes(panel2)
plothypmean(stateratios(ismember(stateratios.SessionIdx,ctrlidx),:),1);
xlim(xlim1);
axes(panel2)
plothypmean(stateratios(ismember(stateratios.SessionIdx,rolidx),:),2);
xlim(xlim1);
legendStr=string(injlist);
legend(legendObjects,legendStr);
pairs={[1 2]};
positions={[.9 1]};

ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures/f5');
for ipair=1:numel(pairs)
    pair=pairs{ipair};
    stats=statistics.fieldtrip.TwoConditionComparison( ...
        ctdatas{pair(1)},ctdatas{pair(2)});
    try
        stat=stats.getClusterBasedIndepSamplesT;
        axes(panel2);
        pos=positions{ipair};
        stat.plot(pos);
        % ff.save(sprintf('Running-Theta-Freq-%s-%s-%s',join(string(states_sub)), ...
        %     join(string(blockstr)),join(string(pair))));
    catch ME
    end
end

ff.save(sprintf('Running-Theta-Power-%s-%s',join(string(states_sub)), ...
            join(string(blockstr))));
%%
function ax1 = plothypmean(tbl,order)
    ax = gca;
    ax1 = axes(Position = ax.Position);
    ax1.Position(2) = ax1.Position(2) + ax1.Position(4)+(order-1)*ax1.Position(4) / 10;
    ax1.Position(4) = ax1.Position(4) / 10;
    vars = ["A-WAKE", "Q-WAKE", "SWS", "REM"];
    for v = vars
        tbl.(v{1})(isnan(tbl.(v{1}))) = 0;
    end
    tbl.awakeRatio = (tbl.("A-WAKE") + tbl.("Q-WAKE")) ./ ...
        (tbl.("A-WAKE") + tbl.("Q-WAKE") + tbl.SWS + tbl.REM);
    a = pivot(tbl, Columns = "Session", Rows = "ZTCenter", ...
        DataVariable = "awakeRatio", IncludeMissingGroups = false);
    data = table2array(a(:, 2:end))';
    data(data == 0) = nan;
    data = smoothdata(data, 2, 'gaussian', 10, 'omitmissing');
    imagesc(hours(a.ZTCenter), 1, mean(data, "omitmissing"));
    cl = flipud(othercolor('RdBu9', 20));
    colormap(gca, cl)
    clim([0 1]);
    ax1.XAxis.Visible = "off";
    ax1.YAxis.Visible = "off";
    ax1.XTickLabel; drawnow;
end