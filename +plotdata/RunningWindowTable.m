classdef RunningWindowTable

    properties
        Table
    end

    methods
        function obj = RunningWindowTable(table)
            %RUNNINGWINDOWTABLE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Table = table;
        end
        function objs = getUnique(obj,var)
            mthds = unique(obj.Table.(var));
            obj1=obj;
            for i=1:numel(mthds)
                obj1.Table=obj.Table(obj.Table.(var)==mthds(i),:);
                objs(i)=obj1;
            end
        end
        function obj = getVarVal(obj,var,val)
            obj.Table = obj.Table(obj.Table.(var)==val,:);
        end
        function ch = getChannel(obj,xvar,yvar)
            ch.t = obj.Table.(xvar)';
            ch.val=obj.Table.(yvar)';
        end
        function p = plot(obj,xvar,yvar)
            chan= obj.getChannel(xvar,yvar);
            if isduration(chan.t)
                chan.t=hours(chan.t);
            end
            p=plot(chan.t,smoothdata(...
                smoothdata(chan.val,'movmedian',5,'includemissing'),...
                'gaussian',5));
        end
        function t = getTime(obj)
            t= unique(obj.Table.ZTCenter);
        end
        function [mat,time1,tbl] = getMerged(obj,group,val)
            ts= obj.getUnique(group);
            time1=obj.getTime';
            tbl=[];
            mat=nan(numel(ts),numel(time1));
            for it=1:numel(ts)
                t=ts(it);
                t1=unique(t.Table(:,{'Session','Condition','ps_method'}));
                idx=ismember(time1,t.Table.ZTCenter');
                mat(it,idx)=t.Table.(val)';
                tbl=[tbl;t1];
            end
        end
    end
end

