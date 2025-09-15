classdef ErrorBarTable
    %ERRORBARTABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x
        ymatrix
    end
    
    methods
        function obj = ErrorBarTable(x,ytable)
            %ERRORBARTABLE Construct an instance of this class
            %   Detailed explanation goes here
            if isduration(x)
                x=hours(x);
            end
            obj.x = x;
            obj.ymatrix = ytable;
        end
        
        function p = plotShaded(obj)
            p = plotdata.ErrorBarPlotProps( ...
                shadedErrorBar(obj.x,obj.getyMean,obj.getyErr));
            p.mainLine.LineWidth=1.5;
        end
        function n = getyN(obj)
            n = sum(~isnan(obj.ymatrix),1);
        end
        function std1 = getyMean(obj)
            std1 = mean(obj.ymatrix,1,"omitmissing");
        end
        function std1 = getySTD(obj)
            std1 = std(obj.ymatrix,1,"omitmissing");
        end
        function std1 = getyErr(obj)
            std1 = obj.getySTD./sqrt(obj.getyN);
        end
    end
end

