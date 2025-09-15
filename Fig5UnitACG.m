% Load required data files
load('Scripts/Theta/tblresBouts.mat')
load('Scripts/Theta/resc.mat')
load('Scripts/Theta/stateratios.mat')

% Create session factory and get sessions table
sf=experiment.SessionFactory;
st=sf.getSessionsTable;

% Load cell metrics sessions and figure factory
fuce=experiment.plot.unit.FiguresUnitCE;

%%
% Clear all variables except the ones listed
clearvars -except tblresBouts resc tses sf st fuce tmres stateratios

% Get cell metrics sessions
cms=fuce.CellMetricsSessions;

% Load unit info table
load("Scripts/Theta/Unit/unitinfotable.mat");

% Create figure factory instance for saving figures
ff=logistics.FigureFactory.instance('./Scripts/Theta');

% Loop over all cell metric sessions
for is=1:numel(cms)
    colors=colororder; % Get color order for plots
    cm=cms(is); % Current cell metric session
    ses=find(ismember(st.Filepath, cm.Session.SessionInfo.baseFolder)); % Find session index
    % Get state ratios for current session
    stateses=stateratios(ismember(stateratios.Session,sprintf('ses%d', ...
        ses(is))),:);
    % Find bouts for current session
    boutidx=ismember(tblresBouts.Session,sprintf('ses%d',ses));
    acgres=[];
    if any(boutidx)
        % Find resc rows for current session
        boutidx2=ismember(resc.Session,sprintf('ses%d',ses));
        tblrescses=resc(boutidx2,:);
        % Get units for current session
        unitses=clusterInfoTable(clusterInfoTable.session==ses,:);
        for iunit=1:height(unitses)
            unit=unitses.id(iunit); % Current unit id
            % Get ACG and theta modulation table for session and unit
            [acgall,tmall]=getACGforSession(ses,unit,tblrescses);
            if ~isempty(acgall)
                % Concatenate results if not empty
                if exist('acgres','var')&&~isempty(acgres)
                    try
                        acgres=[acgres; acgall];
                    catch ME
                        % Ignore errors in concatenation
                    end
                else
                    acgres=acgall;
                end
                clf % Clear figure
                tl1=tiledlayout(5,2); % Create tiled layout
                nexttile(1,[4 1])
                ax=plotacgtime(acgall,tmall); % Plot ACG time
                hypax=plothyp(stateses); % Plot hypnogram
                nexttile(tl1,9,[1 1]);
                % Plot theta frequency and modulation
                plot(hours(tmall.ZTCenter),tmall.lmfreq,Color=colors(1,:));
                hold on
                freq2=(4*tmall.w)./(5-4*tmall.d);
                plot(hours(tmall.ZTCenter),freq2,Color=colors(2,:));
                plot(hours(tmall.ZTCenter),tmall.w,'-',Color=colors(2,:));
                axdetail=gca;
                yyaxis(axdetail,"right");
                plot(hours(tmall.ZTCenter),tmall.a./tmall.b,Color='k');
                ylabel('theta modulation')
                axdetail.YLim=[0 .1];
                linkaxes([ax hypax axdetail],'x'); % Link x axes
                ff.save(sprintf('acgtime')) % Save figure
            end
        end
    end
    % Save ACG results for all units in session
    fname=sprintf("Scripts/Theta/Unit/acgtable-ses%d-allunits.mat",ses);
    save(fname,"acgres")
end

% Function to plot ACG time heatmap and theta modulation
function [ax]=plotacgtime(acgall,tmall)
    normalizewind=[50 200]/1000; % Normalization window in seconds
    acg=pivot(acgall,"DataVariable","count",Columns="ZTCenter",Rows="t"); % Pivot table
    mat0=table2array(acg(:,2:end)); % Get matrix of counts
    t=acg.t; % Time vector
    varszt=acg.Properties.VariableNames; % Variable names (ZTCenter)
    for iv=2:numel(varszt)    
        str=varszt{iv};
        x(iv-1)=str2double(str(1:end-3)); % Extract ZTCenter as numeric
    end
    mats1=smoothdata(mat0,1,"gaussian",50,"omitmissing"); % Smooth data along time
    mats1n=normalize(mats1,1,"zscore"); % Normalize across time
    mats2n=smoothdata(mats1n,2,"gaussian",5,"omitmissing"); % Smooth across ZTCenter
    imagesc(x,t,mats2n); % Plot heatmap
    idx=t>min(normalizewind)&t<max(normalizewind); % Index for normalization window
    ax=gca;
    ax.CLim=[min(mats2n(idx,:),[],"all") max(mats2n(idx,:),[],"all")]; % Set color limits
    ax.XLim=[0 5];
    ax.YLim=[0 .35];
    colormap(ax,"parula");
    hold on
    % Calculate and plot theta modulation as scatter
    thetamodulation=normalize(tmall.a./tmall.b,"range",[2 5]);
    tfreq=x;
    colors=colororder;
    freq=tmall.lmfreq;
    scatter(tfreq,1./freq,thetamodulation,repmat(colors(1,:),length(tfreq),1),'filled',AlphaData=.5)
    scatter(tfreq,2./freq,thetamodulation,repmat(colors(1,:),length(tfreq),1),'filled',AlphaData=.5)
    freq2=(4*tmall.w)./(5-4*tmall.d);
    scatter(tfreq,1./freq2,thetamodulation,repmat(colors(2,:),length(tfreq),1),'filled',AlphaData=.5)
    scatter(tfreq,2./freq2,thetamodulation,repmat(colors(2,:),length(tfreq),1),'filled',AlphaData=.5)
end

% Function to get ACG and theta modulation table for a session and unit
function [acgall,tmall]=getACGforSession(ses, unit, tblrescses)
    acgall=[];
    tmall=[];
    for iwind=1:height(tblrescses)
        fnameacg=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d-acg.mat", ...
                    ses,iwind); % Filename for ACG
        fnametm=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d.mat", ...
            ses,iwind); % Filename for theta modulation
        if isfile(fnameacg)
            s=load(fnameacg);
            acg=s.acg;
            s=load(fnametm);
            tm=s.tm;
        end
        idx=ismember(acg.Properties.VariableNames,string(unit)); % Find unit column
        if any(idx)
            acg1=acg(:,idx); % Extract unit ACG
            tm1=tm(tm.id==unit,:); % Extract unit theta modulation
            acg1.Properties.VariableNames={'count'};
            acg1.ZTCenter=repmat(tblrescses.ZTCenter(iwind),height(acg1),1);
            acg1.Session=repmat(ses,height(acg1),1);
            acg1.unit=repmat(unit,height(acg1),1);
            acg1.t=acg.t;
            if exist('acgall','var')&&~isempty(acgall)
                try
                    acgall=[acgall; acg1];
                    tmall=[tmall;tm1];
                catch ME
                    % Ignore errors in concatenation
                end
            else
                acgall=acg1;
                tmall=tm1;
            end       
        end
    end
end

% Function to plot hypnogram (wake/sleep ratio) above main plot
function ax1=plothyp(tbl)
    ax=gca;
    ax1=axes(Position=ax.Position);
    ax1.Position(2)=ax1.Position(2)+ax1.Position(4);
    ax1.Position(4)=ax1.Position(4)/10;
    awakeSum=sum([tbl.("A-WAKE") tbl.("Q-WAKE")],2,"omitmissing");
    sleepSum=sum([tbl.("SWS") tbl.("REM")],2,"omitmissing");
    awakeSum(isnan(awakeSum))=0;sleepSum(isnan(sleepSum))=0;
    awakeRatio=awakeSum./(awakeSum+sleepSum); % Calculate awake ratio
    imagesc(hours(tbl.ZTCenter),1,awakeRatio');
    cl=flipud(othercolor('RdBu9',200));
    colormap(gca,cl)
    clim([0 1]);
    ax1.XAxis.Visible="off";
    ax1.YAxis.Visible="off";
    ax1.XTickLabel;drawnow;
end