clear all % Clear workspace

% Create session factory for aged experiment
sf=experiment.SessionFactoryAged;
st=sf.getSessionsTable; % Get sessions table

f1=figure(1);clf; % Create and clear figure
tiledlayout('vertical'); % Use vertical tile layout

% Define what to plot and plotting order
plots=categorical({'cf','power','bw','Offset','Knee','Exponent','r_squared','error'});
whattoplot=plots(1); % Select first plot type
plotorder=[2 9 3 14 5 1 10 4 13 6 7 12 15 8 11 16]; % Order of sessions to plot

% FOOOF settings
settings.peak_width_limits=[2 10];
settings.max_n_peaks=8;
settings.aperiodic_mode='knee';%'knee''fixed'
f_range=[1 40]; % Frequency range
dur=minutes(30);slide=minutes(3); % Window duration and slide step

foooftable_sum=[]; % Initialize summary table

% Loop over sessions in plotorder
for i=plotorder
    nexttile % Move to next tile in layout
    ses1=sf.getSessions(i); % Get session object
    ff=logistics.FigureFactory.instance(fullfile(ses1.getBasePath, 'cache/')); % Figure factory
    dlfp=ses1.getDataLFP; % Get LFP data
    sdd=dlfp.getStateDetectionData; % Get state detection data
    cache1=cache.Manager.instance(fullfile(ses1.getBasePath, 'cache/cacheconf.mat')); % Cache manager
    
    % Try to load FOOOF table from cache
    foooftable=cache1.get('FooofTable');
    if isempty(foooftable)
        % If not cached, compute FOOOF table
        thc=cache1.get('ThetaChannel');
        if isempty(thc)
            thc=dlfp.getChannelTimeDataHard.getChannel(sdd.getThetaChannelID);
            [cache1,key]=cache1.hold(thc,'ThetaChannel');
        end
        ss=sdd.getStateSeries;
        aw=ss.getState('AWAKE'); % Get awake state
        awc=thc.getTimeWindow(aw); % Get theta channel during awake
        ticd=awc.TimeIntervalCombined;
        zt=ticd.getZeitgeberTime; % Zeitgeber time
        st1=hours(round(hours((ticd.getStartTimeZT))*20)/20); % Start time rounded
        en1=hours(round(hours((ticd.getEndTimeZT))*20)/20); % End time rounded
        starts=st1:slide:(en1-dur+slide); % Sliding window starts
        ends=starts+dur; % Sliding window ends
        slidingWindows=array2table([starts' (starts'+dur/2) ends'], ...
            VariableNames={'Start','Center','End'}); % Table of windows
        
        % Initialize tables for results
        Peaks=array2table(nan(height(slidingWindows),4), ...
            VariableNames={'cf','power','bw','powerAbs'});
        Aperiodic=array2table(nan(height(slidingWindows),3), ...
            VariableNames={'Offset', 'Knee','Exponent'});
        Fit1=array2table(nan(height(slidingWindows),2), ...
            VariableNames={'r_squared', 'error'});
        
        % Loop over sliding windows
        for iwin=1:height(slidingWindows)
            winabs=[slidingWindows(iwin,:).Start slidingWindows(iwin,:).End]+zt; % Absolute window times
            try
                awcwindow=awc.getTimeWindow(winabs); % Get windowed data
                fooof=awcwindow.getPSpectrumWelch.getFooof(settings,f_range); % Compute FOOOF
                pk=struct2table(fooof.getPeak([5 10])); % Get peak parameters
                Peaks(iwin,:)=pk(1,:);
                Aperiodic(iwin,:)=array2table(fooof.fooof_results.aperiodic_params);
                Fit1(iwin,:)=array2table([fooof.fooof_results.r_squared fooof.fooof_results.error]);
                f2=figure(Visible="off");clf;fooof.plot;fooof.plotPeaks;ax=gca;ax.YLim=[3 5.5];
                text(0.1,0.1,['ZT:' num2str(hours(slidingWindows(iwin,:).Center))],'Units','normalized');
                ff.save('fooofs')
                close(f2);
            catch ME
                % Handle specific errors, ignore others
                if strcmp(ME.identifier, ...
                        'signal:welchparse:invalidSegmentLength')||...
                        strcmp(ME.identifier, ...
                        'MATLAB:structRefFromNonStruct')
                else
                    throw(ME)
                end
            end
        end
        foooftable=[slidingWindows Peaks Aperiodic Fit1]; % Combine results
        [cache1,key]=cache1.hold(foooftable,'FooofTable'); % Save to cache
        clear cache1;
    end
    % Optionally filter by r_squared
    % idx=foooftable.r_squared<.98;
    % foooftable(idx,4:9)=array2table(nan(sum(idx),6));
    
    figure(f1);
    ax(i)=gca; % Get current axes
    x=hours(foooftable.Center); % X axis: window centers
    y=foooftable.(char(whattoplot)); % Y axis: selected metric
    ym=medfilt1(y,dur/slide,'includenan'); % Median filter
    yms=smoothdata(ym,'gaussian',dur/slide,'includemissing'); % Gaussian smoothing
    colors=colororder;
    p=plot(x,y); % Plot raw data
    p.LineWidth=.2;
    p=plot(x,yms); % Plot smoothed data
    p.LineWidth=2;
    p.Color=colors(1,:);
    
    % Set axis limits based on metric
    switch whattoplot
        case 'cf'
            ax(i).YLim=[5.5 8.5];
        case 'power'
            ax(i).YLim=[.1 .8];
        case 'bw'
            ax(i).YLim=[1.5 5];
        case 'Offset'
            ax(i).YLim=[6 8.5];
        case 'Knee'
            ax(i).YLim=[3 5.5];
        case 'Exponent'
            ax(i).YLim=[2 3];
        case 'r_squared'
            ax(i).YLim=[3 5.5];
        case 'error'
            ax(i).YLim=[3 5.5];
    end
    ax(i).XLim=[-3 21]; % X axis limits
    
    % Plot state series and block info
    ss=sdd.getStateSeries;
    ss_zt=ss.getZTCorrected;
    ss_zt.plot([0 .1]+.1);
    text(0.02,1-.02,ses1.toStringShort,'Units','normalized',VerticalAlignment='top');
    bl=ses1.getBlock;bl.plot(ax(i),[0 .1]-.02);

    % Add session number to table
    sesno=i;
    session_tbl=array2table(repmat(sesno,height(foooftable),1),VariableNames={'Session'});
    foooftable_sum=[foooftable_sum; [session_tbl foooftable]]; % Append to summary
end

linkaxes(ax,'xy'); % Link axes for zooming

% Save figures in different sizes
ff=logistics.FigureFactory.instance('./Scripts/Aged/Theta');
f=gcf;f.Position(3)=200;drawnow;getkeyval
ff.save([char(whattoplot) ...
    'Narrrow']);
f=gcf;f.Position(3)=800;
ff.save([char(whattoplot) ...
    'Wide'])

% Save summary table
save("Scripts/Aged/Theta/fooofRunning.mat","foooftable_sum");

%%
clear;
f=figure(1);clf;tiledlayout("vertical");
f.Units="centimeters";f.Position(3)=15;f.Position(4)=19;

% Create three panels for plotting
pnl1=nexttile;hold on;grid on; 
pnl2=nexttile;hold on;grid on; 
pnl3=nexttile;hold on;grid on; 

sf=experiment.SessionFactoryAged;
colors=colororder;

% Load summary table
s=load("Scripts/Aged/Theta/fooofRunning.mat");foooftable_sum=s.foooftable_sum;
sessions=unique(foooftable_sum.Session);
st=sf.getSessionsTable;
st.Condition=categorical(st.Condition);
st.Injection=categorical(st.Injection);
st.animal=categorical(st.animal);

% Filter out sessions with certain injections
noinj=~ismember(st.Injection,categorical({'ROF','SAL'}));
conditions=unique(st.Condition);
noinjwind=ismember(foooftable_sum.Session,st.SessionNo(noinj));
foooftable_sum=foooftable_sum(noinjwind,:);

% Select block to analyze
block=categorical("PRE");

clear dataMat;
% Loop over conditions
for icond=1:numel(conditions)
    cond=conditions(icond);
    sesidx=st.Condition==cond&~ismember(st.Injection,categorical({'ROF','SAL'}));
    ses_sub=st.SessionNo(sesidx);
    ses_cond=[];
    % Loop over sessions in condition
    for ises=1:numel(ses_sub)
        sesno=ses_sub(ises);
        ses1=sf.getSessions(sesno);
        bl_time=ses1.getBlockZT(block);
        ses_tbl=foooftable_sum(foooftable_sum.Session==sesno,:);
        bl_win_idx=ses_tbl.Start>bl_time(1)&ses_tbl.End<bl_time(2);
        ses_tbl_win=ses_tbl(bl_win_idx,:);
        cf1=ses_tbl_win.cf;
        ses_cond=[ses_cond;ses_tbl_win];
    end
    t=unique(ses_cond.Center); % Unique window centers
    res_mat=nan(numel(ses_sub),numel(t)); % Initialize result matrix
    for ises=1:numel(ses_sub)
        sesno=ses_sub(ises);
        ses_tbl=ses_cond(ses_cond.Session==sesno,:);
        idx=ismember(t,ses_tbl.Center);
        res_mat(ises,idx)=ses_tbl.cf;
    end
    dat=data.basic.ChannelTime([], t', res_mat); % Create ChannelTime object
    matp2=dat.getMedianFiltered(minutes(30)); % Median filter
    dataf{icond}=matp2.getGaussianFiltered(minutes(30)); % Gaussian filter
    axes(pnl1);
    p=dataf{icond}.plotChannels(colors(icond,:)); % Plot channels
    axes(pnl2);
    p=dataf{icond}.plotErrorBar(colors(icond,:)); % Plot error bars

    % Prepare time vector for sleep fraction
    t=unique(ses_cond.Center);t=hours(t);
    t=min(t)-.5:.05:max(t)+.5;
    t=round(t*20)/20;
    res_mat=nan(numel(ses_sub),numel(t));
    for ises=1:numel(ses_sub)
        sesno=ses_sub(ises);
        ses1=sf.getSessions(sesno);
        ss=ses1.getStateSeries;
        bl_time=ses1.getBlockZT(block);
        ss1=ss.getWindow(time.ZT(bl_time));
        sr=ss1.getStateRatios(minutes(30),minutes(3));
        frac=sr.getSleepFractionPreceeding;
        frac.ZTtime=round(hours(frac.ZTtime)*20)/20;
        idx=ismember(t,frac.ZTtime');
        try
            res_mat(ises,idx)=frac.sleepFraction;
        catch ME
        end
    end
    dat=data.basic.ChannelTime([], hours(t), res_mat);
    matp2=dat.getMedianFiltered(minutes(30));
    datas{icond}=matp2.getGaussianFiltered(minutes(30));
    axes(pnl3);
    p=datas{icond}.plotErrorBar(colors(icond,:)); % Plot sleep fraction error bars
end

axes(pnl1);
linkaxes([pnl1 pnl2],'xy'); % Link axes for zooming
pnl1.DataAspectRatio
ylim(pnl1,[5.6 8.4])
linkaxes([pnl1 pnl3],'x');
ylim(pnl3,[0 1])
xlim([min(t) max(t)]);
pnl1.DataAspectRatio=[1 .5 1];
pnl2.DataAspectRatio=[1 .5 1];
pnl3.DataAspectRatio=[diff(pnl1.YLim) .5 1];
xticks([pnl1,pnl2,pnl3],-3:25)
axes(pnl2);
axs=plothypmean(datas); % Plot hypnogram mean
linkaxes([pnl2 axs(1) axs(2)],'x');

% Statistical comparison between conditions
st=statistics.fieldtrip.TwoConditionComparison(dataf{1},dataf{2});
stat=st.getClusterBasedIndepSamplesT;
axes(pnl2);
stat.plot()
axes(pnl1);
stat.plot()

st=statistics.fieldtrip.TwoConditionComparison(datas{1},datas{2});
stat=st.getClusterBasedIndepSamplesT;
axes(pnl3);
stat.plot()

% Save figure
ff=logistics.FigureFactory.instance('./Scripts/Aged/Theta');
ff.save(sprintf('Runningwindow-%s.mat',block))

% Function to plot hypnogram mean
function axs = plothypmean(chans)
    ax = gca;
    for ichan = 1:numel(chans)
        ax1 = axes(Position = ax.Position);
        ax1.Position(2) = ax1.Position(2) + ax1.Position(4)+ax1.Position(4) / 10*(ichan-1);
        ax1.Position(4) = ax1.Position(4) / 10;
        data=chans{ichan};
        data1=1-data.Data;
        zt=hours(data.Time);
        imagesc(zt, 1, mean(data1, "omitmissing"));
        cl = flipud(othercolor('RdBu9', 20));
        colormap(gca, cl)
        clim([0 1]);
        ax1.XAxis.Visible = "off";
        ax1.YAxis.Visible = "off";
        ax1.XTickLabel; drawnow;
        axs(ichan)=ax1;
    end
end
