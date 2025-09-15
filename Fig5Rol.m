%% get SD NSD REMs
sf=experiment.SessionFactoryJ;
st=sf.getSessionsTable;
conditions=categorical({'NSD','SD'});
% blocks=categorical({'PRE'});states=categorical({'REM','AWAKE','QWAKE'});
% blocks=categorical({'NSD'});states=categorical({'REM','AWAKE','QWAKE'});
% blocks=categorical({'TRACK'});states=categorical({'AWAKE','QWAKE'});
selected_ses=1:17;idx_ses=ismember(st.ID,selected_ses);
youwantsleep=unique(st.SLEEP);idx_sleep=ismember(st.SLEEP,youwantsleep);
youwantinjection=unique(st.INJECTION);idx_inj=ismember(st.INJECTION,youwantinjection);
youwantanimal=unique(st.ANIMAL);idx_animal=ismember(st.ANIMAL,youwantanimal);
idx=idx_animal&idx_inj&idx_sleep&idx_ses;
st_sub=st(idx,:);
%%
blocks=categorical({'PRE','RUN','NSD','RS'});
blockyouwant=blocks([3 4 2 1]);
states=categorical({'AWAKE','QWAKE','REM'});
stateyouwant=states(1:2);
%%
tblres=[];
for ises= 1:height(st_sub)
    st_ses=st_sub(ises,:);
    for iblock=1:numel(blockyouwant)
        blockwin=blockyouwant(iblock);

        for istate=1:numel(stateyouwant)
            state=stateyouwant(istate);

            fname=sprintf(['Scripts/Theta/matfiles-rol/tblres(' ...
                'bl%s-st%s-ses%d).mat'],blockwin, ...
                state,st_ses.ID);
            if isfile(fname)
                load(fname)
            else
                ses=sf.getSessions(st_ses.ID);
                win_sd=ses.getBlock(blockwin)+ses.SessionInfo.Date;
                if numel(win_sd)>1
                    sest=experiment.SessionTheta(ses);
                    ss=sest.getStates.getWindow(win_sd);
                    [st_time,st_tbl]=ss.getState(state);
                    height1=height(st_tbl);
                    if height1>0
                        thetaratio=sest.getThetaRatioBuzcode;
                        sw=sest.getSlowWaveBroadbandBuzcode;
                        try
                            emg=sest.getEMG.getTimeWindow(win_sd);
                        catch ME
                        end
                        % try
                        %     spd=sest.getSpeed.getTimeWindow(win_sd);
                        % catch ME
                        % end
                        th_sd=sest.getThetaChannel.getTimeWindow(win_sd);

                        clear sesstr
                        sesstr.Condition=repmat(st_ses.SLEEP,height1,1);
                        sesstr.Block=repmat(blockwin,height1,1);
                        sesstr.Session=repmat(st_ses.ID,height1,1);
                        sesstr.State=repmat(state,height1,1);
                        sesstr.Speed=neuro.basic.Channel.empty(height1,0);
                        sesstr.Signal=neuro.basic.ChannelTheta.empty(0,1);
                        sesstr.EMG=neuro.basic.ChannelsThreshold.empty(0,1);
                        sesstr.ThetaRatio=neuro.basic.ChannelsThreshold.empty(0,1);
                        sesstr.SWSlope=neuro.basic.ChannelsThreshold.empty(0,1);
                        for ibout=1:height(st_tbl)
                            twin=[st_tbl(ibout,:).AbsStart st_tbl(ibout,:).AbsEnd];
                            sesstr.Signal(ibout,1)=th_sd.getTimeWindow(twin);
                            sesstr.EMG(ibout,1)=emg.getTimeWindow(twin);
                            % try
                            %     sesstr.Speed(ibout,1)=spd.getTimeWindow(twin);
                            % catch ME
                            %     sesstr.Speed(ibout,1)=[];
                            % end
                            sesstr.ThetaRatio(ibout,1)=thetaratio.getTimeWindow(twin);
                            sesstr.SWSlope(ibout,1)=sw.getTimeWindow(twin);
                            sesstr.Signal(ibout,1)=th_sd.getTimeWindow(twin);
                            % ch.getFrequencyThetaInstantaneous
                        end
                        tbl_ses=[struct2table(sesstr) st_tbl(:,{'AbsStart','AbsEnd'})];
                        % plottbl(tbl_ses,ses)
                        save(fname,'tbl_ses')
                    end
                end
            end
            tblres=[tblres; tbl_ses];
        end
    end
end
save(sprintf('Scripts/Theta/matfiles-rol/tblres-%s-%s.mat', ...
    join(string(stateyouwant)),join(string(blockyouwant))),'tblres','-v7.3');
%% Get Fooofs for every bout and save them in a table
clearvars -except tblres stateyouwant blockyouwant sf
s=load(sprintf('Scripts/Theta/matfiles-rol/tblres-%s-%s.mat', ...
    join(string(stateyouwant)),join(string(blockyouwant))));
tblres=s.tblres;
durs=tblres.AbsEnd - tblres.AbsStart;
idx=durs>seconds(3);
tblres(~idx,:)=[];
tbl_theta_peak=[];
%%
for ibout=1:height(tblres)
    bout=tblres(ibout,:);
    signal1=bout.Signal;
    try
        ps=signal1.getPSpectrumWelch;
        fooof1=ps.getFooof;
        pk=fooof1.getPeak([5 10]);
        bout1=bout(:,{'AbsStart','AbsEnd','Block','Condition','State','Session'});
        bout=[bout1 struct2table(pk)];
        tbl_theta_peak=[tbl_theta_peak; bout];
    catch ME

    end
end
%% add ZT times
sess=sf.getSessions(1:17);
zts=datetime.empty(0,1);
for ibout=1:height(tbl_theta_peak)
    bout=tbl_theta_peak(ibout,:);
    ses=sess(bout.Session);
    zts(ibout,1)=ses.SessionInfo.Date+ses.SessionInfo.ZeitgeberTime;
end
ztstart=tbl_theta_peak.AbsStart-zts;
ztend=tbl_theta_peak.AbsEnd-zts;
tbl_theta_peak=[tbl_theta_peak array2table(zts,VariableNames={'ZT'})];
tbl_theta_peak=[tbl_theta_peak array2table(ztstart,VariableNames={'ZTstart'})];
tbl_theta_peak=[tbl_theta_peak array2table(ztend,VariableNames={'ZTend'})];
save(sprintf('Scripts/Theta/matfiles-rol/tbl_theta_peak-%s-%s.mat', ...
    join(string(stateyouwant)),join(string(blockyouwant))),'tbl_theta_peak','-v7.3')
%% Get Running window theta
s=load(sprintf('Scripts/Theta/matfiles-rol/tblres-%s-%s.mat', ...
    join(string(stateyouwant)),join(string(blockyouwant))));tblres=s.tblres;clear s;
%%
blocks=categorical({'PRE','NSD','RUN','RS'});
sf=experiment.SessionFactoryJ;
selected_ses=1:17;
%             selected_ses=[1 2 11 12 21 22 ];
tses=sf.getSessions(selected_ses);
for iblock=1:numel(blocks)
    stateratios=[];
    blockstr=blocks(iblock);
    for ises=1:numel(tses)
        ses=tses(ises);
        % zt=ses.SessionInfo.ZeitgeberTime+ses.SessionInfo.Date;
        blockwin=ses.getBlockZT(blockstr);

        lfp=ses.getDataLFP;
        % zt=ses.getZeitgeberTime;
        % lfp.TimeIntervalCombined=lfp.TimeIntervalCombined.setZeitgeberTime(zt);
        % lfp.TimeIntervalCombined.saveTable;
        sdd=lfp.getStateDetectionData;
        ss=sdd.getStateSeries;
        ss1=ss.getWindow(blockwin+ss.TimeIntervalCombined.getZeitgeberTime);
        rats=ss1.getStateRatios(minutes(30),minutes(3),blockwin);
        a=array2table(repmat(ises,height(rats),1),VariableNames={'SessionIdx'});
        b=array2table(repmat(selected_ses(ises),height(rats),1), ...
            VariableNames={'Session'});
        c=array2table(repmat(categorical(string(ses.SessionInfo.Condition)), ...
            height(rats),1), ...
            VariableNames={'Condition'});
        subres=[a b c rats];
        stateratios=[stateratios;subres];
    end
    save(sprintf('Scripts/Theta/matfiles-rol/stateratios-%s.mat',blockstr), ...
        'stateratios','tses')
end
%% Get running window of fooof theta and save it in a table
blockyouwant=categorical({'PRE','RUN','RS'});%,'NSD'
sf=experiment.SessionFactoryJ;
for iblock=1:numel(blockyouwant)
    blockstr=blockyouwant(iblock);
    s=load(sprintf("Scripts/Theta/matfiles-rol/stateratios-%s.mat",blockstr));
    s_ratios=s.stateratios;
    sessions=unique(s_ratios.Session);
    % states=categorical("REM");
    states_sub=categorical({'AWAKE','QWAKE'});
    %%
    runningWindowTheta=[];
    for ises=1:numel(sessions)
        sesno=sessions(ises);
        tblres_ses=tblres(tblres.Session==sesno,:);
        ses=sf.getSessions(sesno);
        zt=ses.SessionInfo.ZeitgeberTime+ses.SessionInfo.Date;
        ztcenter=tblres_ses.AbsStart+(tblres_ses.AbsEnd-tblres_ses.AbsStart)/2-zt;
        idx_ses=s_ratios.Session==sesno;
        s_sub=s_ratios(idx_ses,:);
        peaktbl=[];
        for isub=1:height(s_sub)
            win1=s_sub(isub,:);
            win2=[win1.ZTStart win1.ZTEnd];
            tblres_ses_win=tblres_ses(ztcenter>win2(1)&ztcenter<win2(2),:);
            tblres_ses_win1=tblres_ses_win(ismember(tblres_ses_win.State, ...
                states_sub),:);
            if height(tblres_ses_win1)>0
                for ibout=1:height(tblres_ses_win1)
                    sig=tblres_ses_win1(ibout,:).Signal;
                    if ibout==1
                        sig_combined=sig;
                    else
                        sig_combined=sig_combined+sig;
                    end
                end
                if sig_combined.getLength>seconds(3)
                    p_welch=sig_combined.getPSpectrumWelch;
                    p_welch_fooof=p_welch.getFooof;
                    pk=p_welch_fooof.getPeak([5 10]);
                else
                    pk.bw=nan;pk.cf=nan;pk.power=nan;
                end
            else
                pk.bw=nan;pk.cf=nan;pk.power=nan;
            end
            peaktbl=[peaktbl; struct2table(pk)];
        end
        st=array2table(repmat(categorical(join(string(states_sub))), ...
            height(s_sub),1), "VariableNames",{'States'});
        tbl_ses=[s_sub(:,{'Session','Condition','ZTStart','ZTCenter', ...
            'ZTEnd'}) peaktbl st];
        runningWindowTheta=[runningWindowTheta; tbl_ses];
    end
    save(sprintf('Scripts/Theta/matfiles-rol/runningWindowTheta-%s-%s.mat', ...
        join(string(states_sub)),join(string(blockstr))), ...
        'runningWindowTheta','-v7.3')
end

%% plot running Windows
blockstr=categorical("PRE");xlim1=[-4 0];
% blockstr=categorical("RUN");xlim1=[0 1];
% blockstr=categorical("NSD");xlim1=[0 5];
% blockstr=categorical("RS");xlim1=[5 10];

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
            cf=tbl_ses.cf;
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
        panel1.DataAspectRatio=[2 1 1];
        panel2.DataAspectRatio=[2 1 1];
        panel1.XTick=0:5;
        panel2.XTick=0:5;
        panel2.YLim=[5.6 8.4];
        linkaxes([panel1 panel2],'xy')
        if blockstr=="RUN"
            xlim1=hours([0 max(x)]);
        else
            xlim1=round(hours([min(x) max(x)]));
        end
        xlim(xlim1);
    end
end
xlabel('ZT (h)');ylabel('Frequency (Hz)')
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
        ff.save(sprintf('Running-Theta-Freq-%s-%s-%s',join(string(states_sub)), ...
            join(string(blockstr)),join(string(pair))));
    catch ME
    end
end

ff.save(sprintf('Running-Theta-Freq-%s-%s',join(string(states_sub)), ...
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