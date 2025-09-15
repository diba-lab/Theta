% Load required data tables
load('Scripts/Theta/tblresBouts.mat')
load('Scripts/Theta/stateratios.mat')

% Create session factory and get sessions table
sf=experiment.SessionFactory;
st=sf.getSessionsTable;

% Create plotting object
fuce=experiment.plot.unit.FiguresUnitCE;

%%
% Define filename for saving/loading processed data
fname="Scripts/Theta/Speed/corr1tbl.mat";
alltbl=[];

% Check if processed data exists, otherwise process and save
if isfile(fname)
    % Load previously processed data
    s=load(fname);T1=s.T1;corr1tbl2=s.corr1tbl2;alltbl=s.alltbl;clear s;
else
    % Filter bouts with high theta power and frequency below 10 Hz
    tblresBouts1=tblresBouts(tblresBouts.PowerFooof>.75,:);
    tblresBouts1=tblresBouts1(tblresBouts1.CentralFrequenyFooof<10,:);
    
    % Create figure factory for saving figures
    ff=logistics.FigureFactory.instance('./Scripts/Theta/Speed');
    sf=experiment.SessionFactory;
    
    % Select session indices to analyze
    selected_ses=[1:2 4:12 14:15 18:23 ];
    tses=sf.getSessions(selected_ses);
    
    % Define colors for plotting
    colors=othercolor('PRGn8',2);
    colorqwake=colors(2,:);colorawake=colors(1,:);
    isplot=true; % Flag for plotting
    
    % Initialize variables for storing results
    tableBoutSpeed=[];
    spdsesbouts=[];
    statessesbouts=[];
    sesidsesbouts=[];
    corr1tbl.ses=[];
    corr1tbl.cond=[];
    corr1tbl.state=[];
    corr1tbl.window=[];
    corr1tbl.p=[];
    corr1tbl.R=[];
    corr1tbl.coef1=[];
    corr1tbl.coef0=[];
    item=1;
    
    % Loop over selected sessions
    for ises=1:numel(tses)
        ses=tses(ises);
        pos=ses.getPosition; % Get position data
        
        if ~isempty(pos)
            tw=ses.getBlock('SD'); % Get block for sleep deprivation
            possd=pos.getWindow(tw); % Get position data for block
            
            if height(possd.data)>1
                % Plot position and speed if enabled
                if isplot
                    fig1=figure(1);clf;fig1.Position=[2562 0 1439/3 1805/3];
                    tl1=tiledlayout(3,3);
                    ax1=nexttile(1,[1 3]);
                    possd.plot
                end
                ax1.YLabel.String='Postition (cm)';
                
                % Get and filter speed data
                spd=possd.getSpeed;
                spd=spd.getMedianFiltered(1).getGaussianFiltered(3);
                
                if isplot
                    ax2= nexttile(4,[1 3]);
                    spd.plot
                    linkaxes([ax1 ax2],'x')
                end
                ax2.XLim=[0 5];
                
                % Select bouts for current session
                subtbl=tblresBouts1( ...
                    strcmp(tblresBouts1.Session,sprintf('ses%d',selected_ses(ises))) ...
                    ,:);
                speed=nan(height(subtbl),1);
                et=spd.TimeIntervalCombined.getEndTime;
                if isplot
                    axes(ax2);
                end
                ax2.YLabel.String='Speed (cm/s)';
                st=spd.TimeIntervalCombined.getStartTime;
                
                % Loop over bouts and calculate mean speed for each
                for itbl=1:height(subtbl)
                    win=[subtbl(itbl,:).startAbs subtbl(itbl,:).endAbs];
                    if(st<win(1)&&et>win(1))||(st<win(2)&&et>win(2))
                        ab=spd.getTimeWindow(win);
                        speed(itbl)=mean(ab.Values);
                    end
                    % Highlight bout windows on plot
                    if isplot
                        x=hours(win-spd.TimeIntervalCombined.getZeitgeberTime);
                        y=[ax2.YLim(1) ax2.YLim(1) ax2.YLim(2) ax2.YLim(2)];
                        coloridx=find(ismember({'QWAKE','AWAKE'},subtbl(itbl,:).state));
                        patch([x fliplr(x)],y,colors(coloridx,:),FaceAlpha=.2,EdgeColor='none')
                    end
                end
                
                % Add speed to bout table
                speedtbl=array2table(speed,VariableNames={'Speed'});
                subtbl=[subtbl speedtbl];
                
                % Remove rows with NaN values
                subtbl=subtbl(~any(isnan([subtbl.Speed, ...
                    subtbl.CentralFrequenyFooof])'),:);
                tableBoutSpeed=[tableBoutSpeed; subtbl];
                
                % Store speed windows and states for further analysis
                win=[subtbl.startAbs subtbl.endAbs];
                spdsesbouts=[spdsesbouts;spd.getTimeWindow(win)];
                stidx=subtbl.state=="QWAKE";
                win=[subtbl(stidx,:).startAbs subtbl(stidx,:).endAbs];
                spdsesbouts=[spdsesbouts; spd.getTimeWindow(win)];
                stidx=subtbl.state=="AWAKE";
                win=[subtbl(stidx,:).startAbs subtbl(stidx,:).endAbs];
                spdsesbouts=[spdsesbouts; spd.getTimeWindow(win)];
                statessesbouts=[statessesbouts;"WAKE"; "QWAKE";"AWAKE"];
                sesidsesbouts=[sesidsesbouts; repmat(selected_ses(ises),3,1)];
                
                % Define time windows for analysis
                twindows={hours([0 5/3]) hours([5/3 10/3]) hours([10/3 5])};
                st=hours(hours(subtbl.startAbs-ses.SessionInfo.ZeitgeberTime...
                    -ses.SessionInfo.Date));
                en=hours(hours(subtbl.endAbs-ses.SessionInfo.ZeitgeberTime...
                    -ses.SessionInfo.Date));
                subtbl=[subtbl array2table(st,VariableNames="StartZT")...
                    array2table(en,VariableNames="EndZT")];
                
                % For each time window, calculate correlation and plot
                if isplot
                    for iwin=1:numel(twindows)
                        twin=twindows{iwin};
                        idx=st>twin(1)&en<twin(2);
                        if sum(idx)>8
                            tbltoplot=subtbl(idx,:);
                            [R,P] = corrcoef(tbltoplot.Speed,tbltoplot.CentralFrequenyFooof);
                            % Perform linear regression
                            p = polyfit(tbltoplot.Speed,tbltoplot.CentralFrequenyFooof, 1);
                            % Fit a first-degree polynomial (line) to the data
                            x_fit = linspace(min(tbltoplot.Speed), max(tbltoplot.Speed), ...
                                numel(tbltoplot.Speed));
                            y_fit = polyval(p, x_fit); % Generate y-values from the fitted line
                            
                            % Plot the data points and the linear regression line
                            ax3(iwin)=nexttile(6+iwin);hold on;
                            states=sort(unique(subtbl.state),'descend');
                            sizes=seconds(tbltoplot.endAbs-tbltoplot.startAbs);
                            for istate=1:numel(states)
                                stateidx=tbltoplot.state==states(istate);
                                s(istate)=scatter(tbltoplot(stateidx,:).Speed, ...
                                    tbltoplot(stateidx,:).CentralFrequenyFooof,sizes(stateidx), ...
                                    'filled', ...
                                    MarkerFaceAlpha=.7);
                            end
                            hold on;
                            ax3(iwin).XLim=[0 10];ax3(iwin).YLim=[5 10];
                            % Define the different marker sizes
                            sizes = [20, 50, 100];
                            % Create a custom legend for marker sizes
                            legendSizes = {'20 s', '50 s', '100 s'};
                            for i = 1:numel(sizes)
                                scatter(7, max(ax3(iwin).YLim) - (numel(sizes) + i)*0.2, ...
                                    sizes(i), 'filled', 'MarkerFaceColor', 'k', ...
                                    MarkerFaceAlpha=.5);
                                text(7 + 0.4, max(ax3(iwin).YLim) - (numel(sizes) + i)*0.2, ...
                                    legendSizes{i});
                            end
                            % Plot regression line
                            plot(x_fit, y_fit,Color=[.3 .3 .3],LineWidth= 2);
                            % Display regression equation and correlation
                            text(1,6,sprintf('Freq=Speed*%.2f+%.2f, R=%.3f, p=%f', ...
                                p(1),p(2),R(2,1),P(2,1)),"FontWeight","bold", ...
                                "FontSize", 9, "Color",[.3 .3 .3]);
                            xlabel('Speed (cm/s)');ylabel('Frequency (Hz)')
                            legend(s,{'QWAKE','AWAKE'});
                            hold off;
                            % Annotate plot with session info
                            text(0.01,1-0.01,strcat(ses.toStringShort, ...
                                sprintf(' %.2f-%.2f',hours(twin(1)),hours(twin(2)))), ...
                                Units="normalized", VerticalAlignment="top",FontSize=5);
                            % Store correlation results
                            corr1tbl.ses=[corr1tbl.ses;selected_ses(ises)];
                            corr1tbl.cond=[corr1tbl.cond;string(ses.SessionInfo.Condition)];
                            corr1tbl.state=[corr1tbl.state;"WAKE"];
                            corr1tbl.window=[corr1tbl.window;twin];
                            corr1tbl.p=[corr1tbl.p;P(2,1)];
                            corr1tbl.R=[corr1tbl.R;R(2,1)];
                            corr1tbl.coef1=[corr1tbl.coef1;p(1)];
                            corr1tbl.coef0=[corr1tbl.coef0;p(2)];
                        end
                    end
                    % Annotate position and speed plots
                    text(ax1,0.01, 1-0.01,ses.toStringShort,Units='normalized',VerticalAlignment='top',FontSize=5)
                    text(ax2,0.01, 1-0.01,ses.toStringShort,Units='normalized',VerticalAlignment='top',FontSize=5)
                    ff.save(strcat('Fig4-Speed-Freq',ses.getBaseName))
                end
            end
        end
        % Append results for all sessions
        alltbl=[alltbl;subtbl];
    end
    % Convert columns to categorical for easier analysis
    alltbl.Session=categorical(alltbl.Session);
    alltbl.Block=categorical(alltbl.Block);
    alltbl.Condition=categorical(alltbl.Condition);
    % Convert correlation results to table
    corr1tbl2=struct2table(corr1tbl);
    corr1tbl2.cond=categorical(corr1tbl2.cond);
    corr1tbl2.state=categorical(corr1tbl2.state);
    % Create summary table for speed and states
    T1=[array2table(sesidsesbouts,VariableNames={'SessionID'})...
        array2table(categorical(statessesbouts),VariableNames={'States'})...
        array2table(spdsesbouts,VariableNames={'Speed'})...
        ];
    % Save processed data
    save(fname,'T1',"corr1tbl2","alltbl");
end

%%
% Save summary table for further analysis
t1=alltbl(:,{'Condition','Session','state','StartZT','EndZT', ...
    'CentralFrequenyFooof','PowerFooof','Speed'});
save('Scripts/Theta/Speed/theta_speed.mat', 't1')

%%
% Run hierarchical bootstrap analysis (function not shown)
hierarchicalbootstrap

%%
% Plot correlation coefficients and regression results
figure(9);clf;tiledlayout('flow');colors=colororder;
conds=unique(corr1tbl2.cond);
vars={'coef0','coef1','R','p'};
for icond=1:numel(conds)
    cond=conds(icond);
    corrsub=corr1tbl2(corr1tbl2.cond==cond,:);
    for ivar=1:numel(vars)
        if icond==1
            axs(ivar)=nexttile;
        else
            axes(axs(ivar));
        end
        hold on
        var=vars{ivar};
        sess=unique(corrsub.ses);
        winlist=unique(corrsub.window,"rows");
        vals=nan(numel(sess),size(winlist,1));
        for ises=1:numel(sess)
            winplot=corrsub(corrsub.ses==sess(ises),:);
            xcenters=hours(winplot.window(:,1)/5*3)+1;
            x=xcenters+icond*.2-.3;
            vals(ises,xcenters)=winplot.(var);
            p=plot(x,winplot.(var),Color=colors(icond,:));
            p.Marker='.';
            p.MarkerSize=10;
            p.LineWidth=.2;
        end
        xe=(1:size(vals,2))-.3+.2*icond;
        xlim([.5 3.5]);
        title(vars{ivar})
        ye=mean(vals,"omitmissing");
        erre=std(vals,"omitmissing")./sqrt(sum(~isnan(vals)));
        eb=errorbar(xe,ye,erre);
        eb.Color=colors(icond,:);
        eb.Marker='.';
        eb.MarkerSize=20;
        eb.LineWidth=2;
    end
end
ff.save('corr fits')

%% add speed over resc
% Load rescue data and calculate mean speed in 30-min windows
load('Scripts/Theta/resc.mat')
sessions=unique(resc.Session);
windowlength=minutes(30);
tspeed=[];
for ises=1:numel(sessions)
    sessionCode=sessions{ises};
    sessionID=str2double(sessionCode(4:end));
    if ismember(sessionID, unique(T1.SessionID))
        session=sf.getSessions(sessionID);
        rescsub=resc(ismember(resc.Session,sessionCode),:);
        zt=session.SessionInfo.Date+session.SessionInfo.ZeitgeberTime;
        T1sub=T1(T1.SessionID==sessionID,:);
        spdses=nan([height(rescsub)*3 1]);
        stateses=[];
        for iwin=1:height(rescsub)
            window=rescsub(iwin,:);
            ztc=window.ZTCenter+zt;
            windowTime=[ztc-windowlength/2 ztc+windowlength/2];
            for istate=1:height(T1sub)
                spd1=T1sub(istate,:).Speed;
                spd=spd1.getTimeWindow(windowTime);
                spdses((iwin-1)*3+istate)=mean(spd.Values,"includemissing");
            end
            stateses=[stateses; T1sub.States];
        end
        sessionIDs=repmat(sessionID,height(rescsub)*3,1);
        t2=[array2table(sessionIDs,VariableNames={'SessionID'})...
            array2table(stateses,VariableNames={'States'})...
            array2table(spdses,VariableNames={'Speed'})...
            array2table(repelem(rescsub.ZTCenter,3),VariableNames={'ZTCenter'})
            ];
        tspeed=[tspeed;t2];
    end
end
save('Scripts/Theta/tspeed.mat','tspeed')

%% plot tspeed
% Plot mean speed for each session and state over time
load('Scripts/Theta/tspeed.mat')
sf=experiment.SessionFactory;
st=sf.getSessionsTable;
ids=unique(tspeed.SessionID);
states=unique(tspeed.States);
conds=st.Condition(ids);
condlist=unique(conds);
figure;t1=tiledlayout(6,2);
plotcount=1;axs=[];
for icond=1:numel(condlist)
    cond=condlist(icond);
    condids=ids(ismember(conds,cond));
    for icondids=1:numel(condids)
        axs(plotcount)=nexttile((icondids-1)*2+icond);plotcount=plotcount+1;
        sesid=condids(icondids);
        ses=sf.getSessions(sesid);
        dlfp=ses.getDataLFP;
        sdd=dlfp.getStateDetectionData;
        ss=sdd.getStateSeries;
        sestspeed=tspeed(tspeed.SessionID==sesid,:);
        for ist=1:numel(states)
            state_tspeed=sestspeed(ismember(sestspeed.States,states(ist)),:);
            plot(state_tspeed.ZTCenter, state_tspeed.Speed);
            hold on;
        end
        text(0.01,0.01,ses.toStringShort, Units="normalized" ...
            , VerticalAlignment="bottom");
        ax=gca;ax.XLim=hours([0 5]);ax.YLim=[0 7];
        ss.plot([.9 1])
    end
end
linkaxes(axs,'xy');
