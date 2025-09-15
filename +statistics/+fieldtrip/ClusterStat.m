classdef ClusterStat
    %CLUSTERSTAT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Struct
        ClusterTable
    end

    methods
        function obj = ClusterStat(statstruct)
            %CLUSTERSTAT Construct an instance of this class
            %   Detailed explanation goes here
            t=hours(hours(seconds(statstruct.time)));
            int=seconds(unique(diff(statstruct.time)))/2;
            signs={'pos','neg'};
            tbl=[];
            for isign=1:numel(signs)
                sign=signs{isign};
                try
                    clu=struct2table(statstruct.(sprintf('%sclusters',sign)));

                    for iclu=1:height(clu)
                        label=statstruct.(sprintf('%sclusterslabelmat',sign));
                        timeStart(iclu,1)=min(t(label==iclu))-int;
                        timeEnd(iclu,1)=max(t(label==iclu))+int;
                    end
                    if height(clu)>0
                        tblsub=[clu array2table(timeStart) array2table(timeEnd)];
                        if isempty(tbl)
                            tbl=tblsub;
                        else
                            tbl=[tbl;tblsub];
                        end
                    end
                catch ME
                end
            end
            obj.Struct = statstruct;
            obj.ClusterTable=tbl;
        end

        function [] = plot(obj,plotrangeYNorm)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            ax=gca;
            if ~exist('plotrangeYNorm','var')
                plotrangeYNorm=[.9 1];
            end
            plotrangeYabs=ax.YLim(1)+[diff(ax.YLim)*plotrangeYNorm(1)...
                diff(ax.YLim)*plotrangeYNorm(2)];
            tbl=obj.ClusterTable;
            colorRange=[1 .05 .01 .001];
            alpha=[.1 .5 .75 1];
            colors=colororder;
            hold on
            obj.Struct.stat=fillmissing(obj.Struct.stat,'constant',0);
            df=numel(obj.Struct.cfg.design)-2;
            p_vals=[.001 .01 .05];
            t_vals1 = tinv(1-p_vals/2, df);
            t_vals=sort([-t_vals1 t_vals1]);
            % other=othercolor('RdBu9');
            other=othercolor('PiYG9');
            other=othercolor('BrBG9');
            other=othercolor('PuOr9');
            cmap_t=linspace(min(t_vals),max(t_vals),size(other,1));
            other1=other;
            other1(abs(cmap_t)<t_vals1(3),:)=1;
            t_stats=obj.Struct.stat;
            % categorizeFunction = @(x) sum(x > t_vals) + 1;
            % t_stats1 = arrayfun(categorizeFunction, cmap_t);
            y=[plotrangeYabs(2)-diff(plotrangeYabs)/4 plotrangeYabs(2)];
            imagesc(ax,hours(seconds([min(obj.Struct.time) max(obj.Struct.time)])), ...
                y,t_stats,[min(t_vals) max(t_vals)]);
            colormap(ax,other1);
            ax.CLim=[floor(min(t_vals)) ceil(max(t_vals))];
            % cb=colorbar;
            % cb.Ticks=round(t_vals,2);cb.Label.String='t-value';
            for icluster=1:height(tbl)
                cluster=tbl(icluster,:);
                x=hours([cluster.timeStart cluster.timeStart ...
                    cluster.timeEnd cluster.timeEnd]);
                y=[plotrangeYabs(1) plotrangeYabs(1)+diff(plotrangeYabs)/2 ...
                    plotrangeYabs(1)+diff(plotrangeYabs)/2 plotrangeYabs(1)];
                if cluster.clusterstat>0
                    color=other1(end,:);
                else
                    color=other1(1,:);
                end
                alphaval=alpha(find(cluster.prob<colorRange,1,'last'));
                p=patch(x,y,color,'FaceColor','none',"LineWidth",alphaval);
                p.EdgeColor=color;
                text(mean(x),mean(y),sprintf('p=%.3f',cluster.prob), ...
                    "HorizontalAlignment","center",Rotation=45,FontSize=6)
            end
        end
        function report(obj)
            tbl = obj.ClusterTable;
            if isempty(tbl) || height(tbl) == 0
                disp('No significant clusters found.');
                return;
            end

            fprintf('--- Cluster-based Permutation Test Report ---\n');
            for i = 1:height(tbl)
                cluster = tbl(i,:);

                % Determine polarity and related label matrix
                if cluster.clusterstat > 0
                    signLabel = 'Positive';
                    labelmat = obj.Struct.posclusterslabelmat;
                else
                    signLabel = 'Negative';
                    labelmat = obj.Struct.negclusterslabelmat;
                end

                % Time window (converted from hours to seconds)
                tStart = seconds(cluster.timeStart);
                tEnd   = seconds(cluster.timeEnd);

                % Identify t-values and involved channels
                mask = labelmat == i;
                nChans = sum(any(mask, 2));
                tvals = obj.Struct.stat(mask);
                [~, peakIdx] = max(abs(tvals));
                peakT = tvals(peakIdx);

                % Report cluster details
                fprintf(['Cluster #%d (%s): t=[%.3f %.3f] s, clusterstat=%.3f, ' ...
                    'p=%.4g, nChans=%d, peakT=%.3f\n'], ...
                    i, signLabel, tStart, tEnd, ...
                    cluster.clusterstat, cluster.prob, ...
                    nChans, peakT);
            end
        end


    end
end

