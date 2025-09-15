classdef StatePeaksTableRetro<neuro.state.StatePeaksTable
    %STATEPEAKSTABLERETRO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = StatePeaksTableRetro(obj1)
            %STATEPEAKSTABLERETRO Construct an instance of this class
            %   Detailed explanation goes here
            fnames=fieldnames(obj1);
            for i=1:numel(fnames)
                obj.(fnames{i})=obj1.(fnames{i});
            end
        end
         function s = plotscatterRetroState(obj, color,dur,st)
            % Plots retro scatter plot with the specified color.
            IS_RETRO = true;
            s = obj.plotScatterCommon(color, IS_RETRO,dur,st);
        end
        function [s,p] = plotscatterRetro(obj, color,dur)
            % Plots retro scatter plot with the specified color.
            IS_RETRO = true;
            s = obj.plotScatterCommon(color, IS_RETRO,dur);
            % Fit linear regression model
            nan1=isnan(s.XData)|isnan(s.YData);
            x=s.XData(~nan1);y=s.YData(~nan1);
            coeffs = polyfit(x, y, 1);

            % Calculate fitted values using original x values
            yFit = polyval(coeffs, x);

            % Plotting
            hold on; % Keeps the scatter plot visible
            % Generate points for regression line for plotting
            xFit = linspace(min(x), max(x), 100);
            yFitPlot = polyval(coeffs, xFit);
            plot(xFit, yFitPlot, '-k'); % Regression line in black

            % Calculate R-squared
            yMean = mean(y);
            SSres = sum((y - yFit).^2);
            SStot = sum((y - yMean).^2);
            R2 = 1 - SSres/SStot;
            
            % Add R-squared value to the plot
            textLocation = max(x) - 0.4*range(x); % Adjust location as needed
            text(textLocation, min(y) + .2*range(y), sprintf('R^2 = %.2f', R2));

            hold off; % Release the plot
        end

        function obj=getStateWithCertainHistory(obj,state,within,mindur)
            tbl=obj.Table;
            within=minutes(within);
            varstr=sprintf('%s %dmins',string(state),within);
            idx=tbl.(varstr)>seconds(mindur);
            obj.Table=tbl(idx,:);
        end
    end
end