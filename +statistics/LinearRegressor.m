classdef LinearRegressor
    properties
        Model
    end
    
    methods
        function obj = trainModel(obj, dataTable)
            % Train the linear regression model using the given data table
            if width(dataTable) < 2
                error(['Expected a table with at least 1 dependent variable' ...
                    ' and 1 predictor.']);
            end   
            obj.Model = fitlm(dataTable);
        end
        
        function predictions = predict(obj, predictorTable)
            % Use the trained model to make predictions
            if width(predictorTable) ~= width(obj.Model.Coefficients.Estimate) - 1
                error('Mismatch in number of predictors provided for prediction.');
            end
            
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            
            predictions = predict(obj.Model, predictorTable);
        end
        
        function weights = getWeights(obj)
            % Return the weights (coefficients) of the trained model
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            
            weights = obj.Model.Coefficients.Estimate;
        end
        
        function Rsquared = getRsquared(obj)
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            Rsquared = obj.Model.Rsquared.Ordinary;
        end
        function model = getModelTxt(obj)
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            model = evalc('disp(obj.Model.Formula)');
        end
        
        function adjRsquared = getAdjustedRsquared(obj)
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            adjRsquared = obj.Model.Rsquared.Adjusted;
        end
        
        function RMSE = getRMSE(obj)
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            RMSE = obj.Model.RMSE;
        end
        function plotModel(obj)
            % Plot the model weights (coefficients) using a bar plot and 
            % annotate with p-values
            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end

            coefficients = obj.Model.Coefficients.Estimate(2:end);
            predictorNames = obj.Model.CoefficientNames(2:end);
            pvalues = obj.Model.Coefficients.pValue(2:end);

            % Create the bar plot
            bar(predictorNames,coefficients);
            title('Model Weights (Coefficients) with p-values');
            xlabel('Predictors');
            ylabel('Weight Value');
            xtickangle(45);  % Rotate x-axis labels for better visibility
            grid on;
            % modeltxt=disp(obj.Model);
            % text(1-.05,1-.05,modeltxt.plotEffects, ...
            %     VerticalAlignment="top",HorizontalAlignment="right", ...
            %     Units="normalized");
            % Annotate bars with p-values
            for i = 1:length(coefficients)
                text(i, coefficients(i), sprintf('p=%.4f', pvalues(i)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'FontSize', 8);
            end
        end
        function obj = stepwiseSelection(obj, dataTable)
            % Perform stepwise regression for feature selection

            if width(dataTable) < 2
                error(['Expected a table with at least 1 dependent variable' ...
                    ' and 1 predictor.']);
            end

            % Stepwise regression
            obj.Model = stepwiselm(dataTable, 'Criterion', 'AIC');

            % Display the selected model
            disp(obj.Model);
        end
        function plotPredictors(obj)
            % Plot each predictor against the dependent variable with the 
            % model's linear fit

            if isempty(obj.Model)
                error('Model has not been trained. Use trainModel first.');
            end
            dataTable=obj.Model.Variables;
            % Extract dependent variable and predictors
            Y = dataTable{:, end};
            predictors = dataTable{:, 1:end-1};
            numPredictors = size(predictors, 2);

            % Create a figure
            for i = 1:numPredictors
                ax1=nexttile(i);

                % Scatter plot of predictor vs dependent variable
                scatter(predictors(:, i), Y, 'o');
                hold on;

                % Overlay the linear relationship from the model
                [sortedPredictors, sortIdx] = sort(predictors(:, i));
                predictedY = predict(obj.Model, dataTable(sortIdx, :));
                plot(sortedPredictors, predictedY, 'r-', 'LineWidth', 2);

                % Labels and title
                xlabel(dataTable.Properties.VariableNames{i});
                ylabel(dataTable.Properties.VariableNames{end});
                title(sprintf('%s vs. %s', ...
                    dataTable.Properties.VariableNames{end}, ...
                    dataTable.Properties.VariableNames{i}));
                grid on;
                hold off;
                ax1.DataAspectRatio=[.2 1 1];
                ax1.XLim=[0 1];
                ax1.YLim=[5.1 8.4];
                ax(i)=ax1;
            end
            linkaxes(ax,'xy');
        end


        % ... Additional functions for other metrics ...
    end
end
