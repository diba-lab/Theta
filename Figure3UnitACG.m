% fua=experiment.plot.unit.FiguresUnitACG;
% fua.plotACG
load('Scripts/Theta/tblresBouts.mat')
load('Scripts/Theta/resc.mat')
load('Scripts/Theta/stateratios.mat')
sf=experiment.SessionFactory;
st=sf.getSessionsTable;
fuce=experiment.plot.unit.FiguresUnitCE;
%%
clearvars -except tblresBouts resc tses sf st fuce tmres stateratios
cms=fuce.CellMetricsSessions;
load("Scripts/Theta/Unit/unitinfotable.mat");
ff=logistics.FigureFactory.instance('./Scripts/Theta');
for is=1:numel(cms)
    colors=colororder;
    cm=cms(is);
    ses=find(ismember(st.Filepath, cm.Session.SessionInfo.baseFolder));
    stateses=stateratios(ismember(stateratios.Session,sprintf('ses%d', ...
        ses(is))),:);
    boutidx=ismember(tblresBouts.Session,sprintf('ses%d',ses));
    acgres=[];
    if any(boutidx)
        boutidx2=ismember(resc.Session,sprintf('ses%d',ses));
        tblrescses=resc(boutidx2,:);
        unitses=clusterInfoTable(clusterInfoTable.session==ses,:);
        for iunit=1:height(unitses)
            unit=unitses.id(iunit);
            [acgall,tmall]=getACGforSession(ses,unit,tblrescses);
            if ~isempty(acgall)
                if exist('acgres','var')&&~isempty(acgres)
                    try
                        acgres=[acgres; acgall];
                    catch ME

                    end
                else
                    acgres=acgall;
                end
                clf
                tl1=tiledlayout(5,2);
                nexttile(1,[4 1])
                ax=plotacgtime(acgall,tmall);
                hypax=plothyp(stateses);
                nexttile(tl1,9,[1 1]);
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
                linkaxes([ax hypax axdetail],'x');
                ff.save(sprintf('acgtime'))
            end
        end
    end
    fname=sprintf("Scripts/Theta/Unit/acgtable-ses%d-allunits.mat",ses);
    save(fname,"acgres")
end

function [ax]=plotacgtime(acgall,tmall)
normalizewind=[50 200]/1000;
acg=pivot(acgall,"DataVariable","count",Columns="ZTCenter",Rows="t");
mat0=table2array(acg(:,2:end));
t=acg.t;
varszt=acg.Properties.VariableNames;
for iv=2:numel(varszt)    
    str=varszt{iv};
    x(iv-1)=str2double(str(1:end-3));
end
mats1=smoothdata(mat0,1,"gaussian",50,"omitmissing");
mats1n=normalize(mats1,1,"zscore");
mats2n=smoothdata(mats1n,2,"gaussian",5,"omitmissing");
imagesc(x,t,mats2n);
idx=t>min(normalizewind)&t<max(normalizewind);
ax=gca;
ax.CLim=[min(mats2n(idx,:),[],"all") max(mats2n(idx,:),[],"all")];
ax.XLim=[0 5];
ax.YLim=[0 .35];
colormap(ax,"parula");
hold on
% tmall=[];
% for i=1:size(mats2n,2)
%     acg1=neuro.spike.AutoCorrelogram();
%     acg1.Count=mats1(:,i)';
%     acg1.Time=t';
%     tm=acg1.getThetaModulationsTable;
%     acg1.Count=mats2n(:,i)';
%     [freq(i), power(i)]=acg1.getThetaPeaksLocalmax;
%     if isempty(tmall)
%         tmall=tm;
%     else
%         tmall=[tmall;tm];
%     end
% end
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

function [acgall,tmall]=getACGforSession(ses, unit, tblrescses)
acgall=[];
tmall=[];
for iwind=1:height(tblrescses)
    fnameacg=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d-acg.mat", ...
                ses,iwind);
    fnametm=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d.mat", ...
        ses,iwind);
    if isfile(fnameacg)
        s=load(fnameacg);
        acg=s.acg;
        s=load(fnametm);
        tm=s.tm;
    end
    idx=ismember(acg.Properties.VariableNames,string(unit));
    if any(idx)
        acg1=acg(:,idx);
        tm1=tm(tm.id==unit,:);
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

            end
        else
            acgall=acg1;
            tmall=tm1;
        end       
    end
end
end
function ax1=plothyp(tbl)
ax=gca;
ax1=axes(Position=ax.Position);
ax1.Position(2)=ax1.Position(2)+ax1.Position(4);
ax1.Position(4)=ax1.Position(4)/10;
awakeSum=sum([tbl.("A-WAKE") tbl.("Q-WAKE")],2,"omitmissing");
sleepSum=sum([tbl.("SWS") tbl.("REM")],2,"omitmissing");
awakeSum(isnan(awakeSum))=0;sleepSum(isnan(sleepSum))=0;
awakeRatio=awakeSum./(awakeSum+sleepSum);
imagesc(hours(tbl.ZTCenter),1,awakeRatio');
cl=flipud(othercolor('RdBu9',200));
colormap(gca,cl)
clim([0 1]);
ax1.XAxis.Visible="off";
ax1.YAxis.Visible="off";
ax1.XTickLabel;drawnow;
end