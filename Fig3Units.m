% fua=experiment.plot.unit.FiguresUnitACG;
% fua.plotACG
load('Scripts/Theta/tblresBouts.mat')
load('Scripts/Theta/resc.mat')
load('Scripts/Theta/stateratios.mat')
sf=experiment.SessionFactory;
st=sf.getSessionsTable;
fuce=experiment.plot.unit.FiguresUnitCE;
%%
clearvars -except tblresBouts resc tses sf st fuce

cms=fuce.CellMetricsSessions;
clusterInfoTable=[];
for is=1:numel(cms)
    cm=cms(is);
    sesno=find(ismember(st.Filepath, cm.Session.SessionInfo.baseFolder));
    zt=cm.Session.SessionInfo.Date+cm.Session.SessionInfo.ZeitgeberTime;
    boutidx=ismember(tblresBouts.Session,sprintf('ses%d',sesno));
    if any(boutidx)
        tblses=tblresBouts(boutidx,:);
        sa=cm.getSpikeArray;clear cm
        cit=sa.ClusterInfo;
        cit.session=repmat(sesno,height(cit),1);
        clusterInfoTable=[clusterInfoTable;cit];
        sa1=sa.getTimeInterval([tblses.startAbs tblses.endAbs]);clear sa
        boutidx2=ismember(resc.Session,sprintf('ses%d',sesno));
        tblrescses=resc(boutidx2,:);
        shift=unique(diff(tblrescses.ZTCenter));
        dur=tblrescses(1,:).ZTCenter*2;

        for iwind=1:height(tblrescses)
            fnamet=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d.mat", ...
                sesno,iwind);
            fnameacg=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d-acg.mat", ...
                sesno,iwind);
            % if ~isfile(fname)
                ztc=tblrescses(iwind,:).ZTCenter;
                ztb=ztc-dur;zte=ztc+dur;
                absb=ztb+zt;abse=zte+zt;
                sa2=sa1.getTimeInterval([absb abse]);
                acg1=sa2.getAutoCorrelogram.getACGWithCountBiggerThan(5);
                acg1.Info.session=repmat(sesno,height(acg1.Info),1);
                acg2=acg1.getSmooth(20,2);
                % acg2=acg2.getSmooth(3,1);
                tm=acg2.getThetaModulationsTable;
                [tm.lmfreq,tm.lmpower]=acg2.getThetaPeaksLocalmax;
                tm.ZTCenter=repmat(ztc,height(tm),1);
                tm.id=acg2.Info.id;
                tm.session=repmat(sesno,height(tm),1);
                idxneg=acg1.Time<0;
                acg1.Time(idxneg)=[];
                acg1.Count(:,idxneg)=[];
                t=array2table(acg1.Time,VariableNames={'t'});
                unittime=array2table(acg1.Count',VariableNames=string( ...
                    acg1.Info.id));
                acg=[t unittime];
                save(fnamet,'tm');
                save(fnameacg,'acg');
            % end
        end
    end
end
save("Scripts/Theta/Unit/unitinfotable.mat",'clusterInfoTable')
%%
clearvars -except tblresBouts resc tses sf st fuce tmres stateratios
cms=fuce.CellMetricsSessions;
load("Scripts/Theta/Unit/unitinfotable.mat");

% if ~exist('tmres','var')
clear tmres
for is=1:numel(cms)
    cm=cms(is);
    sesno=find(ismember(st.Filepath, cm.Session.SessionInfo.baseFolder));
    boutidx=ismember(tblresBouts.Session,sprintf('ses%d',sesno));
    if any(boutidx)
        boutidx2=ismember(resc.Session,sprintf('ses%d',sesno));
        tblrescses=resc(boutidx2,:);

        for iwind=1:height(tblrescses)
            fname=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d.mat", ...
                sesno,iwind);
            fnameacg=sprintf("Scripts/Theta/Unit/acgtable-ses%d-win%d-acg.mat", ...
                sesno,iwind);
            if isfile(fname)
                s=load(fname);
                tm=s.tm;
            end
            if isfile(fnameacg)
                s=load(fnameacg);
                acg=s.acg;
                acg1=table2array(acg)';
                acg2=acg1(2:end,:);
                acg3 = num2cell(acg2, 2); % Convert each row of the array into a cell
                yourTable = cell2table(acg3, 'VariableNames', {'acg'}); % Create the table
            end
            tm1=[tm yourTable];
            if exist('tmres','var')&&~isempty(tmres)
                tmres=[tmres; tm1];
            else
                tmres=tm1;
            end
        end
    end
end
% end
%%
ff=logistics.FigureFactory.instance('./Scripts/Theta');

f1=figure(1);clf(f1);f1.Position=[2564 -1010 1059 527];
tiledlayout(1,1)
nexttile
tmres1=tmres(tmres.adjrsquare>.9,:);
tmres2=tmres1(zscore(tmres1.a./tmres1.b)<3,:);
tmres2=tmres;
[~,neuronsidx]=ismember( string(int2str([tmres2.session tmres2.id]))...
    ,string(int2str([clusterInfoTable.session clusterInfoTable.id]))...
     );
gr=categorical(clusterInfoTable(neuronsidx,"group").group);
ct=categorical(clusterInfoTable(neuronsidx,"cellType").cellType);
br=categorical(clusterInfoTable(neuronsidx,"brainRegion").brainRegion);
good=gr=="good";
mua=gr=="mua";
int=ct=="Narrow Interneuron"|ct=="Wide Interneuron";
pyr=ct=="Pyramidal Cell";
ca3=br=="CA3";
ca1=br=="CA1";
idx=(int|pyr)&(ca1)&(good|mua);
% idx=true(size(good));
tmres2=tmres2(idx,:);
% tmres2=tmres2(region,:);
plotMedianFreqPerSession(tmres2,resc,clusterInfoTable,st)
f2=figure(2);clf(f2);f2.Position=[2702 -1005 1059 1462];
tl=tiledlayout(5,2);
plotFreqPerUnit(tl,tmres2,resc,clusterInfoTable,st,stateratios)
figure(f2);ff.save('Unitall');
figure(f1);ff.save('Unitall');

%%
tmsumm=neuro.spike.ThetaModulationsSummaryTable(tmres);
sess=unique(tmsumm.thetaModulationsTable.session);
f=figure('Units','centimeters');clf;f.Position(3:4)=[5 8];
for ises=1:numel(sess)
    sesno=sess(ises);
    ses=sf.getSessions(sesno);
    ss=ses.getStateSeries();
    sr=ss.getStateRatios(slidingWindowSize, ...
        slidingWindowLaps,windowZT);
    tss1=tmsumm.thetaModulationsTable(tmsumm.thetaModulationsTable.session==sesno,:);
    units=unique(tss1.id );
    for iunit=1:numel(units)
        figure(f);clf;
        unitnumber=units(iunit);
        slidingWindowSize=minutes(30);
        slidingWindowLaps=minutes(5);
        windowZT=hours([0 5]);
        tl=tiledlayout(6,1);
        t1=nexttile(tl,1,[1 1]);
        sr.plotAwakeFractionPreceeding
        colormap(t1,flipud(othercolor('RdBu8')));
        clim([0 1]);
        yticklabels('');
        yticks('');
        xlim([0 5]);
        xticklabels('')
        t2=nexttile(tl,2,[4 1]);
        tmsumm.plotUnit(sesno,unitnumber);
        xlabel('');
        xticklabels('');
        colormap(t2,"parula");
        t3=nexttile(tl,6,[1 1]);
        tmsumm.plotFreq(sesno,unitnumber);
        xlim([0 5]);
        ylim([5.5 9]);
        ylabel('Freq. (Hz)');
        xlabel('ZT (h)');
        str=sprintf('ses:%d\nunit:%d',sesno,unitnumber);
        dim=[0 .7 .3 .3];
        annotation('textbox', dim, 'String', str, 'EdgeColor', 'none', ...
            'VerticalAlignment', 'top');
        ff.save(sprintf('acg_theta_every_unit_in_ses_%d',sesno));
    end
end
%%

function []=plotMedianFreqPerSession(tbl,resc,clu,st)
colors=colororder;
hold on
sescond=st(ismember(st.SessionNo,unique(tbl.session)),:);
sescond.Condition=categorical(sescond.Condition);
sescond.animal=categorical(sescond.animal);
conds=unique(sescond.Condition);
tbl.freqcalculated=tbl.w*4./(5-tbl.d);

% pt=pivot(tbl,"Columns","session",Rows="ZTCenter",DataVariable="freqcalculated", ...
%     Method="median");
% pt=pivot(tbl,"Columns","session",Rows="ZTCenter",DataVariable="w", ...
%     Method="median");
pt=pivot(tbl,"Columns","session",Rows="ZTCenter",DataVariable="lmfreq", ...
    Method="median");
t=hours(pt.ZTCenter);
arr1=table2array(pt(:,2:end))';
arr2=smoothdata(arr1,2,"movmedian","omitmissing");
arr3=smoothdata(arr2,2,"gaussian","omitmissing");
for icond=1:numel(conds)
    cond=conds(icond);
    ses=sescond.SessionNo(sescond.Condition==cond);
    vars=pt.Properties.VariableNames;
    idx=ismember(vars(2:end),string(ses));
    arr4=arr3(idx,:);
    plot(t, arr4,Color=[colors(icond,:) .4]);
    y=mean(arr4,"omitmissing");
    err=std(arr4,"omitmissing")/sqrt(size(arr4,1));
    seb=shadedErrorBar(t,y,err,'lineprops',{'Color',colors(icond,:)});
    seb.mainLine.LineWidth=2;
end
hold off
end

function []=plotFreqPerUnit(tile1,tbl,resc,clu,st,stateratios)
colors=colororder;
sescond=st(ismember(st.SessionNo,unique(tbl.session)),:);
sescond.Condition=categorical(sescond.Condition);
sescond.animal=categorical(sescond.animal);
conds=unique(sescond.Condition);
tbl.freqcalculated=tbl.w*4./(5-tbl.d);
tbl.thetamodulation=tbl.a./tbl.b;

% pt=pivot(tbl,"Columns",{'session','id'},Rows="ZTCenter",DataVariable="freqcalculated", ...
%     Method="median");
% pt=pivot(tbl,"Columns",{'session','id'},Rows="ZTCenter",DataVariable="w", ...
%     Method="median");
pt=pivot(tbl,"Columns",{'session','id'},Rows="ZTCenter",DataVariable="lmfreq", ...
    Method="median");
ptmodulation=pivot(tbl,"Columns",{'session','id'},Rows="ZTCenter",DataVariable="thetamodulation", ...
    Method="median");
t=hours(pt.ZTCenter);
for icond=1:numel(conds)
    cond=conds(icond);
    statecond=stateratios(categorical(stateratios.Condition)==cond,:);    
    ses=sescond.SessionNo(sescond.Condition==cond);
    vars=pt.Properties.VariableNames;
    idx=find(ismember(vars,string(ses)));
    for ises=1:numel(idx)
        pos=(ises-1)*2+icond;
        unitax(pos)=nexttile(tile1,pos);
        hold on
        units=table2array(table2array(pt(:,idx(ises))))';
        samplearray=table2array(pt(1,2));
        unitids=str2double(samplearray.Properties.VariableNames);
        modulation=table2array(table2array(ptmodulation(:,idx(ises))))';
        unitididx=find(~all(isnan(units),2));
        unitsidssub=unitids(unitididx);
        clu1=clu(clu.session==ses(ises),:);
        clu2=clu1(unitsidssub,:);
        arr1=units;
        arr2=smoothdata(arr1,2,"movmedian",10,"omitmissing");
        arr3=smoothdata(arr2,2,"gaussian",10,"omitmissing");
        plot(t, arr3,Color=[colors(icond,:) .2]);
        y=mean(arr3,"omitmissing");
        err=std(arr3,"omitmissing")./sqrt(sum(~isnan(arr3)));
        seb=shadedErrorBar(t,y,err,'lineprops',{'Color',colors(icond,:)});
        seb.mainLine.LineWidth=2;
        types={'group','cellType'};
        clu2.group=categorical(clu2.group);
        clu2.cellType=categorical(clu2.cellType);
        clu2.cellType(clu2.cellType=="Narrow Interneuron")="Wide Interneuron";
        [gr, ia, ic]=unique(clu2(:,types));
        str=[];
        for igr=1:height(gr)
            for ity=1:numel(types)
                ty=gr(igr,:).(types{ity});
                str=sprintf('%s %s',str,ty);
            end
            numcelltype=sum(ic==igr);
            str=sprintf('%s: %d\n',str,numcelltype);
        end

        %         colorbivar=othercolor('Paired8',8);
        % clu2.cellType=categorical(clu2.cellType);
        % clu2.cellType(clu2.cellType=="Narrow Interneuron")="Wide Interneuron";
        % [gr, ia, ic]=unique(clu2(:,{'cellType','group'}));
        % for igr=1:height(gr)
        %     idxgr=ic==igr;
        %     arr4=arr3(idxgr,:);
        %     plot(t, arr4,Color=[colorbivar(igr,:) .5]);
        %     y=mean(arr4,1,"omitmissing");
        %     err=std(arr4,0,1,"omitmissing")./sqrt(sum(~isnan(arr4)));
        %     seb=shadedErrorBar(t,y,err,'lineprops',{'Color',colorbivar(igr,:)});
        %     seb.mainLine.LineWidth=2;
        % end
        txt1=text(1,1,str,Units="normalized");txt1.FontSize=5;
        txt1.HorizontalAlignment="right";txt1.VerticalAlignment="top";
        statecondses=statecond(ismember(statecond.Session,sprintf('ses%d', ...
            ses(ises))),:);
        ylim([7 8.5])
        hypax(pos)=plothyp(statecondses);
        hold off
    end

end
linkaxes([unitax hypax],'x');
end

function ax1=plothyp(tbl)
ax=gca;
ax1=axes(Position=ax.Position);
ax1.Position(2)=ax1.Position(2)+ax1.Position(4);
ax1.Position(4)=ax1.Position(4)/10;
awakeSum=sum([tbl.("A-WAKE") tbl.("Q-WAKE")],2,'omitmissing');
sleepSum=sum([tbl.("SWS") tbl.("REM")],2,'omitmissing');
awakeSum(isnan(awakeSum))=0;sleepSum(isnan(sleepSum))=0;
awakeRatio=awakeSum./(awakeSum+sleepSum);
imagesc(hours(tbl.ZTCenter),1,awakeRatio');
cl=flipud(othercolor('RdBu9',20));
colormap(gca,cl)
clim([0 1]);
ax1.XAxis.Visible="off";
ax1.YAxis.Visible="off";
ax1.XTickLabel;drawnow;
end