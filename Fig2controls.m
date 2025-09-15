
%% get SD NSD REMs
sf=experiment.SessionFactory;selected_ses=[1:2 4:12 14:15 18:23 ];
st=sf.getSessionsTable;
conditions=categorical({'NSD','SD'});
% blocks=categorical({'PRE'});states=categorical({'REM','AWAKE','QWAKE'});
% blocks=categorical({'NSD'});states=categorical({'REM','AWAKE','QWAKE'});
% blocks=categorical({'TRACK'});states=categorical({'AWAKE','QWAKE'});
blocks=categorical({'NSD'});states=categorical({'REM','AWAKE','QWAKE'});
%% Table Results every bout
tblres=[];
for icond=1:numel(conditions)
    condition=conditions(icond);
    idx_cond=ismember(st.Condition,condition);
    idx_selected=ismember(st.SessionNo,selected_ses);
    selected_ses_cond=st.SessionNo(idx_cond & idx_selected);
    for iblock=1:numel(blocks)
        block=blocks(iblock);
        for istate=1:numel(states)
            state=states(istate);
            for ises=1:numel(selected_ses_cond)
                sesno=selected_ses_cond(ises);
                fname=sprintf(['Scripts/Theta/matfiles/tblres(' ...
                    'bl%s-st%s-ses%d).mat'],block,state,sesno);
                if isfile(fname)
                    load(fname)
                else
                    ses=sf.getSessions(sesno);
                    win_sd=ses.getBlock(block)+ses.SessionInfo.Date;
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
                        try
                            spd=sest.getSpeed.getTimeWindow(win_sd);
                        catch ME
                        end
                        th_sd=sest.getThetaChannel.getTimeWindow(win_sd);
                        %%
                        clear sesstr
                        sesstr.Condition=repmat(condition,height1,1);
                        sesstr.Block=repmat(block,height1,1);
                        sesstr.StateTemp=repmat(state,height1,1);
                        sesstr.Session=repmat(sesno,height1,1);
                        sesstr.State=repmat(state,height1,1);
                        sesstr.Speed=neuro.basic.Channel.empty(0,1);
                        sesstr.Signal=neuro.basic.ChannelTheta.empty(0,1);
                        sesstr.EMG=neuro.basic.ChannelsThreshold.empty(0,1);
                        sesstr.ThetaRatio=neuro.basic.ChannelsThreshold.empty(0,1);
                        sesstr.SWSlope=neuro.basic.ChannelsThreshold.empty(0,1);
                        for ibout=1:height(st_tbl)
                            twin=[st_tbl(ibout,:).AbsStart st_tbl(ibout,:).AbsEnd];
                            sesstr.Signal(ibout,1)=th_sd.getTimeWindow(twin);
                            sesstr.EMG(ibout,1)=emg.getTimeWindow(twin);
                            try
                                sesstr.Speed(ibout,1)=spd.getTimeWindow(twin);
                            catch ME
                                sesstr.Speed(ibout,1)=[];
                            end
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
                tblres=[tblres; tbl_ses];
            end
        end
    end
end
save(sprintf('Scripts/Theta/matfiles/tblres-%s-%s.mat', ...
    join(string(states)),join(string(blocks))),'tblres','-v7.3');
%%
clearvars -except tblres blocks states
if ~exist('tblres','var')
    s=load(sprintf('Scripts/Theta/matfiles/tblres-%s-%s.mat', ...
        join(string(states)),join(string(blocks))));tblres=s.tblres;clear s;
    tblres=neuro.state.StateTable(tblres,experiment.SessionFactory);
end
%% Get Fooofs for every bout and save them in a table
tblres=tblres.getDurationLongerThan(seconds(3));
tbl_theta_peak=tblres.getFooofEstimatedTable;
save(sprintf('Scripts/Theta/matfiles/tbl_theta_peak-%s-%s.mat', ...
    join(string(states)),join(string(blocks))),'tbl_theta_peak','-v7.3')
%% Get Running window theta
s=load(sprintf('Scripts/Theta/matfiles/tblres-%s-%s.mat', ...
    join(string(states)),join(string(blocks))));tblres=s.tblres;clear s;
tblres=neuro.state.StateTable(tblres,experiment.SessionFactory);

%% Get running window of fooof theta and save it in a table

s=load(sprintf("Scripts/Theta/stateratios-%s.mat",blocks));s_ratios=s.stateratios;

%% Get Running window theta
states_sub=categorical({'AWAKE','QWAKE'});
runningWindowTheta = tblres.getRunningWindows(s_ratios,states_sub);
save(sprintf('Scripts/Theta/matfiles/runningWindowTheta-%s-%s.mat', ...
    join(string(states_sub)),join(string(blocks))),'runningWindowTheta','-v7.3')
%%
% calculate it first
clear
sf=experiment.SessionFactory; selected_ses=[1:2 4:12 14:15 18:23 ];
st=experiment.SessionList(sf.getSessions(selected_ses),selected_ses);
st.saveRippleEventPropertiesTable;
%% plot sample
% t_interest(1,1)=hours(2)+minutes(.3);duration=seconds(5);
% t_interest(1,2)=t_interest(1,1)+duration;
% win_interest=time.ZT(t_interest);
% sesno_interes=1;
% st.plotRippleEventProperties(sesno_interes,win_interest);
%%
clear
sf=experiment.SessionFactory; selected_ses=[1:2 4:12 14:15 18:23 ];
st=experiment.SessionList(sf.getSessions(selected_ses),selected_ses);
rippleEventPropertiesTableCombined=st.getRippleEventPropertiesTableCombined;
%%
% Eliminate weak channels
% Initialization
sf = experiment.SessionFactory;
rippleEventPropertiesTableCombined1 = groupfilter( ...
    rippleEventPropertiesTableCombined, {'Session', 'peak'}, ...
    @filter1, "power_envelope");
riptbl = neuro.ripple.RippleTable(rippleEventPropertiesTableCombined1);
%%
ff=logistics.FigureFactory.instance('./Scripts/Theta/ripple');
% First scenario
f=figure(1); clf; tiledlayout(10, 2);f.Position=[2562 -1112 759 2453];
vars1 = {"frequency_wavelet", "frequency_wavelet_nowhiten", "frequency_count", "frequency_count_large"};
varstr1 = {'wavelet 100-250', 'wavelet no whitening 100-250', 'count 120-200', 'count large 100-250'};
processAndPlot(sf, riptbl, vars1, varstr1, 'Frequency (Hz)', [150 200],false);
ff.save('RippleFreq')

% Second scenario
f=figure(2); clf; tiledlayout(10, 2);f.Position=[2562 -1112 759 2453];
vars2 = {"power_envelope", "power_envelope_largefilter"};
varstr2 = {'env 120-200', 'env 100-250'};
processAndPlot(sf, riptbl, vars2, varstr2, 'Power (z-score)', [-1.5 1.5],true);
ff.save('RipplePow')

%%
% Second scenario
riptbl = neuro.ripple.RippleTable(rippleEventPropertiesTableCombined);
f=figure(3); clf; tiledlayout(10, 2);f.Position=[2562 -1112 759 2453];
vars2 = {"SWMax"};
varstr2 = {'SW Amplitude'};
processAndPlot(sf, riptbl, vars2, varstr2, 'SW Amp (z-score)', [-1.5 1.5],true);
ff.save('SWAmp')




%% Get Running window ripple
% Get running window of fooof theta and save it in a table
s=load(sprintf("Scripts/Theta/stateratios-%s.mat",blocks));
s_ratios=s.stateratios;
sessions=unique(s_ratios.Session);
runningWindowRipples=[];
for ises=1:numel(sessions)
    sesno=sessions(ises);
    sesfname=sprintf('Scripts/Theta/matfiles-rip/runningWindowRipple-%d-%s.mat', ...
        sesno, join(string(blocks)));
    if ~isfile(sesfname)
        ses=experiment.SessionRipple(sf.getSessions(sesno));
        idx_ses=s_ratios.Session==sesno;
        s_sub=s_ratios(idx_ses,:);
        s=load(sprintf("Scripts/Theta/matfiles-rip/" + ...
            "rippleEventSignalTable-ses%d.mat",sesno));
        rippleEventSignalTable=neuro.ripple.RippleTable( ...
            s.rippleEventSignalTable);clear s;
        runningWindowRipple=rippleEventSignalTable.getRunningWindows(s_sub);
        save(sesfname,'runningWindowRipple','-v7.3')
    else
        s=load(sesfname);runningWindowRipple=s.runningWindowRipple;clear s;
    end
    runningWindowRipples=[runningWindowRipples; runningWindowRipple];
end
save(sprintf('Scripts/Theta/matfiles-rip/runningWindowRipple-%s.mat', ...
    join(string(blocks))),'runningWindowRipples','-v7.3');
%% plot running Windows Ripple
s=load(sprintf('Scripts/Theta/matfiles-rip/runningWindowRipple-%s.mat', ...
    join(string(blocks))));runningWindowRipples=s.runningWindowRipples;clear s;
% rsub=runningWindowRipples(runningWindowRipples.Channel==1,:); clear
% runningWindowRipples;
rsub=runningWindowRipples;clear runningWindowRipples;
%%
isplot=false;
if isplot
    figure(1);clf;tiledlayout("vertical");
    t1=nexttile;t2=nexttile;t3=nexttile;
end
fname1=sprintf('Scripts/Theta/matfiles-rip/runningWindowRipplesCalc-%s.mat', ...
    join(string(blocks)));
if ~isfile(fname1)
    runningWindowRipplesCal=[];
    for iwin=1:height(rsub)
        line1=rsub(iwin,:);
        pw=neuro.power.PowerSpectrumRipple(line1.p_welch);
        pp=neuro.power.PowerSpectrumRipple(line1.p_ps);
        pc=neuro.power.PowerSpectrumRipple(line1.p_ch);
        categories=categorical({'welch','pspectrum','chronux'});
        subt1=line1(:,{'Session','Condition','ZTStart','ZTCenter','ZTEnd'});
        try
            pwf=pw.getFooof;
            pk_w=pwf.getPeak([120 250]);
            runningWindowRipplesCal=[runningWindowRipplesCal ;...
                [subt1 struct2table(pk_w) ...
                array2table(categories(1),VariableNames="ps_method")]];
            ppf=pp.getFooof;
            pk_p=ppf.getPeak([120 250]);
            runningWindowRipplesCal=[runningWindowRipplesCal ;...
                [subt1 struct2table(pk_p) ...
                array2table(categories(2),VariableNames="ps_method")]];
            pcf=pc.getFooof;
            pk_c=pcf.getPeak([120 250]);
            runningWindowRipplesCal=[runningWindowRipplesCal ;...
                [subt1 struct2table(pk_c) ...
                array2table(categories(3),VariableNames="ps_method")]];
            if isplot
                axes(t1)
                pwf.plot
                pwf.plotPeaks
                axes(t2)
                ppf.plot
                ppf.plotPeaks
                axes(t3)
                pcf.plot
                pcf.plotPeaks
            end
        catch ME
            for ic=1:numel(categories)
                runningWindowRipplesCal=[runningWindowRipplesCal ;...
                    [subt1 array2table(nan(1,3),VariableNames={'cf','power','bw'}) ...
                    array2table(categories(ic),VariableNames="ps_method")]];
            end
        end
    end
    save(fname1,'runningWindowRipplesCal','-v7.3')
else
    s=load(fname1);runningWindowRipplesCal=s.runningWindowRipplesCal;clear s;
end

%%
% runningWindowTheta=runningWindowTheta(runningWindowTheta.Condition=="NSD",:);
colors=colororder;
rWRC=plotdata.RunningWindowTable(runningWindowRipplesCal);
t=rWRC.getTime;
ps_methods=rWRC.getUnique("ps_method");
what="cf";mat=[];
for ips_method=1:numel(ps_methods)
    ps_method=ps_methods(ips_method);
    f=figure(1);clf;f.Units="centimeters"; f.Position(3)=10;f.Position(4)=10;
    tl=tiledlayout('vertical');
    panel1=nexttile;grid on;hold on;panel2=nexttile; grid on;hold on;
    conds=ps_method.getUnique("Condition");
    ctdatas=[];
    for icond=1:numel(conds)
        cond=conds(numel(conds)-icond+1);
        [mat{icond},time,tbl]=cond.getMerged("Session",what);
        times{icond}=time;
        matp1=data.basic.ChannelTime([],times{icond},mat{icond});
        matp2=matp1.getMedianFiltered(minutes(60));
        ctdatas{icond}=matp2.getGaussianFiltered(minutes(60));
        axes(panel1);
        p=ctdatas{icond}.plotChannels(colors(icond,:));
        axes(panel2);
        p=ctdatas{icond}.plotErrorBar(colors(icond,:));
    end
    ps_method_str=string(unique(ps_method.Table.ps_method));
    text(.05,1-.05,ps_method_str,VerticalAlignment="top", ...
        HorizontalAlignment="left",Units="normalized");
    linkaxes([panel1 panel2],'xy')
    xlim(round(hours([min(times{icond}) max(times{icond})])));
    panel2.YLim=[150 200];
    xlabel('ZT (h)');ylabel('Frequency (Hz)')

    st=statistics.fieldtrip.TwoConditionComparison(ctdatas{1},ctdatas{2});
    stat=st.getClusterBasedIndepSamplesT;
    axes(panel2);
    stat.plot()    %%
    ff=logistics.FigureFactory.instance('./Scripts/Theta/ripple');
    ff.save(sprintf('Running-Ripple-Freq-%s-%s', ...
        join(string(blocks)),ps_method_str));

end
%% plot Power spectrum time plot.

clearvars -except fname1 blocks
s=load(fname1);runningWindowRipplesCal=s.runningWindowRipplesCal;clear s;
s=load(sprintf('Scripts/Theta/matfiles-rip/runningWindowRipple-%s.mat', ...
    join(string(blocks))));runningWindowRipples=s.runningWindowRipples;clear s;


%% plot running Windows THETA
states_sub=categorical({'AWAKE','QWAKE'});
% states_sub=categorical({'REM'});
blocks=categorical({'NSD'});align=0;xl=[0 5];
% blocks=categorical({'POST'});align=1;xl=[0 3];
% blocks=categorical({'PRE'});align=0;xl=[-3 0];
% blocks=categorical({'TRACK'});align=1;xl=[0 1.5];
s=load(sprintf('Scripts/Theta/matfiles/runningWindowTheta-%s-%s.mat', ...
    join(string(states_sub)),join(string(blocks))));
st=load(sprintf('Scripts/Theta/stateratios-%s.mat',blocks));
srat=st.stateratios;
% Calculate awake ratio and add it to the stateratios table
awakeSum = sum([srat.("A-WAKE") srat.("Q-WAKE")], 2, "omitmissing");
sleepSum = sum([srat.("SWS") srat.("REM")], 2, "omitmissing");
awakeSum(isnan(awakeSum)) = 0;
sleepSum(isnan(sleepSum)) = 0;
srat.awakeRatio = awakeSum ./ (awakeSum + sleepSum);

runningWindowTheta=s.runningWindowTheta;
% runningWindowTheta=runningWindowTheta(runningWindowTheta.Condition=="NSD",:);
f=figure(1);clf;f.Units="centimeters"; f.Position(3)=20;f.Position(4)=17;
tl=tiledlayout('vertical');
panel1=nexttile;grid on;hold on;panel2=nexttile; grid on;hold on;
conditions=sort(unique(runningWindowTheta.Condition),'descend');
sessions=unique(runningWindowTheta.Session);
colors=colororder;
x=unique(runningWindowTheta.ZTCenter);
x1=x-x(1);
if align
    x2=x1(1):hours(1/20):x1(end);
else
    x2=x(1):hours(1/20):x(end);
end

cfs=nan(numel(sessions),numel(x2));
ctdatas=[];
sratses=[];
for icond=1:numel(conditions)
    cond=conditions(icond);
    sratses{icond}=srat(srat.Condition==cond,:);
    runningWindowTheta1=runningWindowTheta( ...
        runningWindowTheta.Condition==cond,:);
    sub_ses=unique(runningWindowTheta1.Session);
    for ises=1:numel(sub_ses)
        sesno=sub_ses(ises);
        ses=sf.getSessions(sesno);bl_zt=ses.getBlockZT(blocks);
        bl_st_zt=hours(round(hours(bl_zt(1))*20)/20);
        tbl_ses=runningWindowTheta1(runningWindowTheta1.Session==sesno,:);
        if align
            cf=interp1(tbl_ses.ZTCenter-bl_st_zt,tbl_ses.cf,x2);
        else
            cf=interp1(tbl_ses.ZTCenter,tbl_ses.cf,x2);
        end
        cfs(ismember(sessions,sesno),:)=cf';
    end
    sub_cfs=cfs(ismember(sessions,sub_ses),:);

    matp1=data.basic.ChannelTime([],x2,sub_cfs);
    matp2=matp1.getMedianFiltered(minutes(30));
    ctdatas{icond}=matp2.getGaussianFiltered(minutes(30));
    axes(panel1);
    p=ctdatas{icond}.plotChannels(colors(icond,:));
    axes(panel2);
    p=ctdatas{icond}.plotErrorBar(colors(icond,:));
end
linkaxes([panel1 panel2],'xy')
xlim(round(hours([min(x2) max(x2)])));
panel2.YLim=[5.6 8.4];
xlabel('ZT (h)');ylabel('Frequency (Hz)')
panel1.DataAspectRatio=[1 .5 1];
panel2.DataAspectRatio=[1 .5 1];

xlim(xl);

st=statistics.fieldtrip.TwoConditionComparison(ctdatas{1},ctdatas{2});
% try
%     stat=st.getClusterBasedIndepSamplesT;
%     stat.report
%     axes(panel2);
%     stat.plot()
% catch ME
% 
% end
if blocks=="NSD"
    str=load('Scripts\Theta\matfiles\trials-AWAKE QWAKE-PRE.mat','ctdatas');
    ctdatasPRE=str.ctdatas;
    st=statistics.fieldtrip.OneConditionComparison(ctdatas{2});
    vals=ctdatas{2};
    vals1=mean(vals.Data(:,1:10),2,"omitmissing");
    vals=ctdatasPRE{2};
    vals1=mean(vals.Data(:,(141:end)),2,"omitmissing");
    % mSD=mean(vals1);
    % stdSD=std(vals1);
    % vals2=ctdatasPRE{1};
    % vals2=mean(vals2.Data(:,(104:end)),2,"omitmissing");
    % mNSD=mean(vals2);
    % stdNSD=std(vals2);
    try
        stat=st.getClusterBaseddepSamplesT(vals1);
        axes(panel2);
        stat.plot()
        stat.report
    catch ME

    end

    plothypmean(sratses{1},x2,xl,0);
    axes(panel2);
    plothypmean(sratses{2},x2,xl,1.1);
    hold on;
    mean_vals = mean(vals1, 'omitnan');
    std_vals = std(vals1, 'omitnan') / sqrt(numel(vals1));
    axes(panel2);
    colors=colororder;thecolor=colors(2,:);
    errorbar(hours(x2(2)), mean_vals, std_vals, 'Color', thecolor, 'LineWidth', 1.5); hold on;
    plot(hours(x2(2)), mean_vals, 'o', 'MarkerFaceColor', thecolor, 'MarkerSize', 6);
end
ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures/f2');
ff.save(sprintf('Running-Theta-Freq-%s-%s-',join(string(states_sub)), ...
    join(string(blocks))));
%%
clearvars -except tbl_theta_peak states blocks
s=load(sprintf('Scripts/Theta/matfiles/tbl_theta_peak-%s-%s.mat', ...
    join(string(states)),join(string(blocks))));
tbl_theta_peak=neuro.state.StatePeaksTable( ...
    s.tbl_theta_peak,experiment.SessionFactory);clear s;
%% Plot bouts in scatter format
state=2:3; conds=categorical({'NSD','SD'});
cond=conds(1);
tbl_theta_peak1=tbl_theta_peak;
ColorOrRd=linspecer(4,'qualitative');
states=unique(tbl_theta_peak.Table.State);
idx_state=ismember(tbl_theta_peak.Table.State,states(state));
idx_cond_nsd=ismember(tbl_theta_peak.Table.Condition,cond);
tbl_theta_peak1.Table=tbl_theta_peak.Table(idx_state&idx_cond_nsd,:);
figure(3);clf;
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 12.5, 8]); % Set figure size
tbl_theta_peak1.plotscatter(ColorOrRd(1,:))
ylim([5.6 8.4])

ff.save('Scatter_time_freq_power_duration')

% tbl_sleep.plotscatter(ColorOrRd(state(1),:))
%% get Retro Duration table
dur=minutes(round(logspace(log10(1),log10(60),7)));
f_thetapeaktable=sprintf('./Scripts/Theta/matfiles/retro-%s-%s.mat', ...
    join(string(cond)),num2str(minutes(dur)));
if ~isfile(f_thetapeaktable)
    tbl_theta_peak2=tbl_theta_peak1.getRetro(dur);
    save(f_thetapeaktable, "tbl_theta_peak2");
else
    s=load(f_thetapeaktable);
    tbl_theta_peak2=s.tbl_theta_peak2;clear s;
end
%% plot all durations sorted
f1=figure(1);clf;f1.Units="centimeters" ;tiledlayout("vertical");
f1.Position(3)=7.5;f1.Position(4)=36;
ax=[];
for idur=4:numel(dur)
    ax1=nexttile;
    dur1= dur(idur);
    tbl_theta_peak2.plotscatterRetro(ColorOrRd(1,:),minutes(dur1));
    % ax1.DataAspectRatio=[.2 .8 1];
    ax1.XLim=[0 1];ax1.YLim=[5.6 8.4];
    ax1.XTick=0:.2:1;
    ax(idur)=ax1;

end
linkaxes(ax,'xy');
ff.save('FITS')
%% plot comparison REM vs SWS wakings
duration=minutes(1);
durmin=seconds(30);
figure(1);clf;tiledlayout(1,4);ax1=nexttile(1,[1 3]);
ax11=nexttile(4,[1 1]);
figure(2);clf;tiledlayout(1,2);ax2=nexttile;ax3=nexttile;
axes(ax1)
tbl_theta_peak2.Table=...
    tbl_theta_peak2.Table(tbl_theta_peak2.Table.ZTstart>hours(2),:);
tbl_sleep=tbl_theta_peak2.getStateWithCertainHistory( ...
    categorical("REM"),duration,durmin);
s1=tbl_sleep.plotscatter(ColorOrRd(2,:));
hold on;
axes(ax11);
tbl_sleep.plotHistogram(ColorOrRd(2,:));view([90 -90]);hold on
axes(ax2);
tbl_sleep.plotscatterRetro(ColorOrRd(2,:),minutes(duration))
tbl_sleep=tbl_theta_peak2.getStateWithCertainHistory( ...
    categorical("SWS"),duration,durmin);
axes(ax1)
s2=tbl_sleep.plotscatter(ColorOrRd(1,:));
legend([s2 s1],{'SWS','REM'})
axes(ax3);
tbl_sleep.plotscatterRetro(ColorOrRd(1,:),minutes(duration))
axes(ax11);xlim(ax1.YLim);
tbl_sleep.plotHistogram(ColorOrRd(1,:));view([90 -90]);hold on
ylabel('pdf')
% linkaxes([ax1, ax11], 'y');
axes(ax1);ylim([5 10]);
axes(ax11);xlim([5 10]);
linkaxes([ax2, ax3], 'xy');
axes(ax2);ylim([5 8.5]);xlim([.5 1]);
ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures/');
figure(1); 
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 12, 6]); % Adjust figure size and position
ff.save(sprintf('Awakening-REMvsSWS1'));

figure(2); 
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, 12, 6]); % Adjust figure size and position
ff.save(sprintf('Awakening-REMvsSWS2'));
%%
dur2=minutes(round(logspace(log10(1),log10(30),4)));
interest=categorical({'SWS','REM'});rest=categorical({'A-WAKE','Q-WAKE'});
%%
tbl_theta_peak2=tbl_theta_peak1.getRetroAlt( ...
    dur2,interest,rest);
save(sprintf('./Scripts/Theta/matfiles/retro-%s-%s-%s-%s.mat', ...
    join(string(cond)),num2str(minutes(dur2)), ...
    join(string(interest)),join(string(rest))), "tbl_theta_peak2");
%%
figure(2),clf,tiledlayout(numel(dur2),numel(interest));
ax=[];
for idur=1:numel(dur2)
    for ist=1:numel(interest)
        ax(idur)=nexttile;
        dur1= dur2(idur);
        st=interest(ist);
        tbl_theta_peak2.plotscatterRetroState( ...
            ColorOrRd(state(1),:),minutes(dur1),st);
    end
end
linkaxes(ax,'xy');
%%
% cf=7;others=[14:17 19:22 24:27 29:32 34:37 39:42 44:47];
cf=7;others=[13:5:43];
istate=1;
% cf=7;others= [14:5:44]+istate-1;
lr=tbl_theta_peak2.getLinearRegressor(cf,others);
lr.getWeights
figure(istate);clf;tiledlayout('flow','TileSpacing','compact');
lr.plotPredictors();
ylim([5 10])
nexttile
lr.Model.plotEffects
modeltxt=sprintf(['Root Mean Squared Error: %.3f\n R-squared: %.3f\n' ...
    'Adjusted R-Squared: %.3f'], ...
    lr.getRMSE,lr.getRsquared,lr.getAdjustedRsquared);
text(1-.05,.05,modeltxt, ...
    VerticalAlignment="bottom",HorizontalAlignment="right", ...
    Units="normalized",FontSize=8);
ff=logistics.FigureFactory.instance('./Scripts/Theta/Figures');
str=join(lr.Model.CoefficientNames(2:end));
ff.save(sprintf('GLM-%s',str{1}))

function []=plottbl(tbl,ses)
for ibout=1:height(tbl)
    tbl1=tbl(ibout,:);
    sig=tbl1.Signal;
    ps=sig.getPSpectrumWelch;
    a=ps.getFooof();
    a.plot;
end
end
function [b]=filter1(x)
b=false(size(x));
b1=find(ismember(x,max(x,[],"includemissing")),1);
b(b1)=true;
if ~any(b)
    b(1)=true;
end
end
% get Track Awake and Qwake
% Function to handle the repetitive plotting

function processAndPlot(sf, riptbl, vars, varstr, ylabelStr, ylimVals,iszscore)
st = sf.getSessionsTable;
sess = riptbl.getSessions;
st = st(ismember(st.SessionNo, sess), :);
conds = sort(unique(st.Condition), "descend");
for ic = 1:numel(conds)
    st1 = st(st.Condition == conds(ic), :);
    sess1 = st1.SessionNo;
    for is = 1:numel(sess1)
        ax(ic, is) = nexttile(is * 2 - (2 - ic));
        sesno = sess1(is);
        ses1=sf.getSessions(sesno);
        win_interest = ses1.getBlockZT("NSD");
        ses = riptbl.getSession(sesno);
        ses = ses.removeArtifacts;
        ses = ses.getTimeWindow(win_interest);
        hold on;
        xlabel('ZT (h)');
        ylabel(ylabelStr);
        for ivar = 1:numel(vars)
            exidx = ismember(ses.Table.(vars{ivar}), ...
                [min(ses.Table.(vars{ivar})) max(ses.Table.(vars{ivar}))]);
            ses.Table(exidx, vars{ivar}) = {nan};
            if iszscore
                tblsub = ses.Table.(vars{ivar});
                tblsub(~isnan(tblsub), 1) = zscore(tblsub(~isnan(tblsub)));
                ses.Table.(vars{ivar}) = tblsub;
            end
            p1 = ses.plotRunningAverage(vars{ivar}, win_interest, ...
                minutes(6), minutes(3));
        end
        txt1=text(0,1.5,ses1.toStringShort);
        ax(ic,is).YGrid="on";
    end
end
legend(varstr);
linkaxes([ax], 'xy');
xlim([0 5]);
ylim(ylimVals);
end

function ax1 = plothypmean(tbl,time,xlim,move)
ax = gca;
ax1 = axes(Position = ax.Position);
ax1.Position(2) = ax1.Position(2) + ax1.Position(4)+move*(ax1.Position(4) / 10);
ax1.Position(4) = ax1.Position(4) / 10;
a = pivot(tbl, Columns = "Session", Rows = "ZTCenter", ...
    DataVariable = "awakeRatio", IncludeMissingGroups = false);
data = table2array(a(:, 2:end))';
data(data == 0) = nan;
data = smoothdata(data, 2, 'gaussian', 10, 'omitmissing');
imagesc(time, 1, mean(data, "omitmissing"));
cl = flipud(othercolor('RdBu9', 20));
colormap(gca, cl)
clim([0 1]);
ax1.XAxis.Visible = "off";
ax1.YAxis.Visible = "off";
ax1.XTickLabel; drawnow;
ax1.XLim=hours(xlim);

end