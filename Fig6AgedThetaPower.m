clear all % Clear workspace

% Create session factory for aged experiments
sf=experiment.SessionFactoryAged;
st=sf.getSessionsTable; % Get sessions table

f1=figure(1);clf; % Create and clear figure
tiledlayout('vertical'); % Use vertical tile layout

% Define what to plot and plot order
plots=categorical({'cf','power','bw','Offset','Knee','Exponent','r_squared','error'});
whattoplot=plots(1); % Select what to plot (first entry)
plotorder=[2 9 3 14 5 1 10 4 13 6 7 12 15 8 11 16]; % Order of sessions to plot

% FOOOF settings
settings.peak_width_limits=[2 10];
settings.max_n_peaks=8;
settings.aperiodic_mode='knee';%'knee''fixed'
f_range=[1 40]; % Frequency range
dur=minutes(30);slide=minutes(3); % Window duration and slide
foooftable_sum=[]; % Initialize summary table

% Loop over sessions in plotorder
for i=plotorder
    nexttile % Move to next tile in layout
    ses1=sf.getSessions(i); % Get session object
    ff=logistics.FigureFactory.instance(fullfile(ses1.getBasePath, 'cache/')); % Figure factory for saving
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
        st1=hours(round(hours((ticd.getStartTimeZT))*20)/20); % Start time (rounded)
        en1=hours(round(hours((ticd.getEndTimeZT))*20)/20); % End time (rounded)
        starts=st1:slide:(en1-dur+slide); % Sliding window starts
        ends=starts+dur; % Sliding window ends
        slidingWindows=array2table([starts' (starts'+dur/2) ends'], ...
            VariableNames={'Start','Center','End'}); % Table of windows
        
        % Initialize tables for FOOOF results
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
                Peaks(iwin,:)=pk(1,:); % Store peak
                Aperiodic(iwin,:)=array2table(fooof.fooof_results.aperiodic_params); % Store aperiodic params
                Fit1(iwin,:)=array2table([fooof.fooof_results.r_squared fooof.fooof_results.error]); % Store fit stats
                f2=figure(Visible="off");clf;fooof.plot;fooof.plotPeaks;ax=gca;ax.YLim=[3 5.5];
                text(0.1,0.1,['ZT:' num2str(hours(slidingWindows(iwin,:).Center))],'Units','normalized');
                ff.save('fooofs') % Save figure
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
    % Plot results
    figure(f1);
    ax(i)=gca;
    x=hours(foooftable.Center); % X axis: window center times
    y=foooftable.(char(whattoplot)); % Y axis: selected parameter
    ym=medfilt1(y,dur/slide,'includenan'); % Median filter
    yms=smoothdata(ym,'gaussian',dur/slide,'includemissing'); % Gaussian smoothing
    colors=colororder;
    p=plot(x,y); % Raw data
    p.LineWidth=.2;
    p=plot(x,yms); % Smoothed data
    p.LineWidth=2;
    p.Color=colors(1,:);
    % Set axis limits based on parameter
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
    ax(i).XLim=[-3 21];
    ss=sdd.getStateSeries;
    ss_zt=ss.getZTCorrected;
    ss_zt.plot([0 .1]+.1); % Plot state series
    text(0.02,1-.02,ses1.toStringShort,'Units','normalized',VerticalAlignment='top'); % Session label
    bl=ses1.getBlock;bl.plot(ax(i),[0 .1]-.02); % Plot block
    
    sesno=i;
    session_tbl=array2table(repmat(sesno,height(foooftable),1),VariableNames={'Session'}); % Session number table
    foooftable_sum=[foooftable_sum; [session_tbl foooftable]]; % Append to summary
end

linkaxes(ax,'xy'); % Link axes for all plots

% Save figures
ff=logistics.FigureFactory.instance('./Scripts/Aged/Theta');
f=gcf;f.Position(3)=200;drawnow;getkeyval
ff.save([char(whattoplot) ...
    'Narrrow']);
f=gcf;f.Position(3)=800;
ff.save([char(whattoplot) ...
    'Wide'])
save("Scripts/Aged/Theta/fooofRunning.mat","foooftable_sum"); % Save summary data

%%
clear;
f=figure(1);clf;tiledlayout("vertical");
f.Units="centimeters";f.Position(3)=15;f.Position(4)=19;
% Create three panels for plotting
pnl1=nexttile;hold on;grid on; pnl2=nexttile;hold on;grid on; pnl3=nexttile;hold on;grid on;
sf=experiment.SessionFactoryAged;
colors=colororder;
s=load("Scripts/Aged/Theta/fooofRunning.mat");foooftable_sum=s.foooftable_sum; % Load summary data
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

% Define block and x-axis limits for plotting
block=categorical("POST");xlim1=[6 11];
block=categorical("TRACK");xlim1=[5 6];
block=categorical("NSD");xlim1=[0 5];
block=categorical("PRE");xlim1=[-3 0];

clear dataMat;
% Loop over conditions
for icond=1:numel(conditions)
    cond=conditions(icond);
    sesidx=st.Condition==cond&~ismember(st.Injection,categorical({'ROF','SAL'}));
    ses_sub=st.SessionNo(sesidx); % Sessions for this condition
    ses_cond=[];
        idx=ismember(t,ses_tbl.Center);
        res_mat(ises,idx)=ses_tbl.power;
    end
    dat=data.basic.ChannelTime([], t', res_mat);
    matp2=dat.getMedianFiltered(minutes(30));
    dataf{icond}=matp2.getGaussianFiltered(minutes(30));
    axes(pnl1);
    p=dataf{icond}.plotChannels(colors(icond,:));
    axes(pnl2);
    p=dataf{icond}.plotErrorBar(colors(icond,:));


    t=unique(ses_cond.Center);t=hours(t);
    t=min(t)-.5:.05:max(t)+.5;
    t=round(t*20)/20;
    %
    % t=hours(round(hours(hours(0):minutes(3):hours(5))*20)/20);
    % t=hours(round(hours(hours(0):minutes(3):hours(5))*20)/20)+hours(6);
    % t=hours(round(hours(hours(0):minutes(3):hours(11))*20)/20)+hours(11);
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
    % axes(pnl1);
    % p=datas{icond}.plotChannels(colors(icond,:));
    axes(pnl3);
    p=datas{icond}.plotErrorBar(colors(icond,:));
end

axes(pnl1);
linkaxes([pnl1 pnl2],'xy');
pnl1.DataAspectRatio
ylim(pnl1,[0.2 1.0])
linkaxes([pnl1 pnl3],'x');
ylim(pnl3,[0 1])
xlim(xlim1);
pnl1.DataAspectRatio=[4 .5 1];
pnl2.DataAspectRatio=[4 .5 1];
pnl3.DataAspectRatio=[diff(pnl1.YLim)*4 .5 1];
xticks([pnl1,pnl2,pnl3],-3:25)
axes(pnl2);
axs=plothypmean(datas);
linkaxes([pnl2 axs(1) axs(2)],'x');
xlim(xlim1);

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


ff=logistics.FigureFactory.instance('./Scripts/Aged/Theta');
ff.save(sprintf('RunningwindowPower-%s.mat',block))

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
