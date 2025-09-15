classdef ErrorBarPlotProps
    %ERRORBARPLOTPROPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mainLine
        patch
        edge
    end
    
    methods
        function obj = ErrorBarPlotProps(props)
            %ERRORBARPLOTPROPS Construct an instance of this class
            %   Detailed explanation goes here
            obj.mainLine = props.mainLine;
            obj.patch = props.patch;
            obj.edge = props.edge;
        end
        
        function obj = setColor(obj,color)
            obj.mainLine.Color=color;
            obj.patch.FaceColor=color;
            for i=1:numel(obj.edge)
                obj.edge(i).Color=color;
            end
        end
    end
end

