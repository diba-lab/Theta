classdef StatePeaksTable
    % Represents a table of state peaks and related operations.

    properties
        Table % Main data table.
        sessionFactory
    end

    methods
        function obj = StatePeaksTable(tbl,sessionFactory)
            % Constructor: Initializes with table 'tbl'.
            if nargin>0
                obj.Table = tbl;
            end
            if nargin>1
                obj.sessionFactory = sessionFactory;
            end
        end

        function obj = getRetro(obj, dur)
            % Extracts retro data based on the provided duration.
            sf = obj.sessionFactory; % Create session factory.
            uniqueSessions = unique(obj.Table.Session);
            resultTbl = table(); % Initialize empty result table.

            for ises = 1:numel(uniqueSessions)
                sessionNumber = uniqueSessions(ises);
                sessionData = obj.Table(obj.Table.Session == sessionNumber,:);
                ses = sf.getSessions(sessionNumber);
                stateSeries = ses.getStateSeries;
                for ibout = 1:height(sessionData)
                    boutData = obj.processBout(sessionData(ibout,:), ...
                        stateSeries, dur);
                    resultTbl = [resultTbl; boutData]; % Append processed data.
                end
            end
            obj.Table = resultTbl; % Assign processed data to main table.
            obj=neuro.state.StatePeaksTableRetro(obj);
        end
        function obj = getRetroAlt(obj, dur,interest,rest)
            % Extracts alternative retro data based on the new parameter.
            % Here, 'someNewParameter' represents any additional data or settings you need.

            sf = obj.sessionFactory; % Create session factory instance.
            uniqueSessions = unique(obj.Table.Session);
            resultTbl = table(); % Initialize empty result table.

            for ises = 1:numel(uniqueSessions)
                sessionNumber = uniqueSessions(ises);
                sessionData = obj.Table(obj.Table.Session == sessionNumber,:);
                ses = sf.getSessions(sessionNumber);
                stateSeries = ses.getStateSeries;

                for ibout = 1:height(sessionData)
                    try
                        % Call the new processing function or modify the
                        % existing one with new logic
                        boutData = obj.processBoutAlt( ...
                            sessionData(ibout,:), stateSeries,dur, ...
                            interest,rest);
                        resultTbl = [resultTbl; boutData]; % Append processed data.
                    catch ME
                        warning(ME.message);
                    end
                end
            end

            obj.Table = resultTbl; % Assign processed data back to main table.
            obj=neuro.state.StatePeaksTableRetro(obj);
        end
        function s = plotscatter(obj, color)
            % Plots scatter plot with the specified color.
            IS_RETRO = false;
            s = obj.plotScatterCommon(color, IS_RETRO);
        end
        function s = plotHistogram(obj, color)
            h=histogram(obj.Table.cf,linspace(5,10,20),Normalization="pdf");
            h.FaceColor=color;
        end
        function lr = getLinearRegressor(obj,var,predictors)
            % Instantiate the LinearRegressor class
            lr = statistics.LinearRegressor();

            % Train the model
            tbl=obj.Table(:,[predictors var]);
            lr = lr.trainModel(tbl);
            % lr = lr.stepwiseSelection(tbl);
        end
        function cvResults = performKFoldCV(obj, dataTable, k)
            % Perform k-fold cross-validation on the model

            if isempty(obj.Model)
                error(['Model has not been trained. Use trainModel or ' ...
                    'stepwiseSelection first.']);
            end

            if width(dataTable) < 2
                error(['Expected a table with at least 1 dependent ' ...
                    'variable and 1 predictor.']);
            end

            % Check if k is greater than 1
            if k <= 1
                error('Number of folds k must be greater than 1.');
            end

            % Prepare the dataset
            responseVarName = dataTable.Properties.VariableNames{1};
            predictors = dataTable(:, 2:end);
            response = dataTable(:, 1);

            % Create a partition for the data
            cvPartition = cvpartition(response.(responseVarName), 'KFold', k);

            % Initialize variables to store results
            testMSEs = zeros(cvPartition.NumTestSets, 1);

            % Perform k-fold cross-validation
            for i = 1:cvPartition.NumTestSets
                % Training/testing indices for this fold
                trainData = dataTable(training(cvPartition, i), :);
                testData = dataTable(test(cvPartition, i), :);

                % Fit the model on the training data
                tempModel = fitlm(trainData);

                % Evaluate the model on the testing data
                yTest = testData{:, 1};
                predictions = predict(tempModel, testData(:, 2:end));

                % Calculate and store the mean squared error
                testMSEs(i) = mean((yTest - predictions).^2);
            end

            % Output the results
            cvResults.averageTestMSE = mean(testMSEs);
            cvResults.testMSEs = testMSEs;
        end
    end

    methods (Access = protected)
        function s = plotScatterCommon(obj, color, isRetro,dur,st)
            % Common method for plotting scatter based on flag isRetro.
            sortedTbl = sortrows(obj.Table, "ZTstart");
            time = hours((sortedTbl.ZTstart + sortedTbl.ZTend) / 2);
            duration = seconds(sortedTbl.ZTend - sortedTbl.ZTstart);
            if isRetro
                if ~exist('st','var')
                    s = scatter(sortedTbl.(sprintf('dur%dmins',dur)), ...
                        sortedTbl.cf, duration, "filled", ...
                        'AlphaData', sortedTbl.power, 'MarkerFaceAlpha', "flat", ...
                        'MarkerFaceColor', color); hold on;
                    xlabel('Sleep Proportion')
                else
                    s = scatter(sortedTbl.(sprintf('dur%dmins%s',dur,st)), ...
                        sortedTbl.cf, duration, "filled", ...
                        'AlphaData', sortedTbl.power, 'MarkerFaceAlpha', "flat", ...
                        'MarkerFaceColor', color); hold on;
                    xlabel(sprintf('Proportion %s',st))
                end
            else
                % Add scatter plot for the main data
                s = scatter(time, sortedTbl.cf, duration, "filled", ...
                    'AlphaData', sortedTbl.power, 'MarkerFaceAlpha', "flat", ...
                    'MarkerFaceColor', color); hold on;

                % Add legend points
                legendTime = [4, 4, 4]-.5;
                legendCf = [8.5, 9, 9.5];
                legendPower = [0.5, 1, 2];
                legendDuration = [50, 50, 50];
                scatter(legendTime, legendCf, legendDuration, "filled", ...
                    'AlphaData', legendPower, 'MarkerFaceAlpha', "flat", ...
                    'MarkerFaceColor', color);

                     % Add legend points
                legendTime = [4, 4, 4];
                legendCf = [8.5, 9, 9.5];
                legendPower = [2, 2, 2];
                legendDuration = [50, 100, 200];
                scatter(legendTime, legendCf, legendDuration, "filled", ...
                    'AlphaData', legendPower, 'MarkerFaceAlpha', "flat", ...
                    'MarkerFaceColor', color);

                % Add legend
                legend('Data Points', 'Legend Points', 'Location', 'best');
                xlim([min(time) max(time)]);
                xlabel('ZT (h)');
            end
            ylabel('Frequency (Hz)');
            grid on;
            vars={'Block','Condition','State'};
            str1="";
            for ivar=1:numel(vars)
                a=join(string(unique(sortedTbl.(vars{ivar}))),'-');
                if ivar==1
                    str1=a;
                else
                    str1=join([str1,a],'|');
                end
            end
            if isRetro
                str1=join([str1,sprintf("%d min",dur)],'|');
            end
            text(.05,1-.05,str1, ...
                Units="normalized", ...
                VerticalAlignment="top", ...
                HorizontalAlignment="left");
           
        end
        function boutData = processBout(~, bout, stateSeries, dur)
            % Processes a bout of data to extract relevant information.
            timeSeriesTbl = table();
            for idur = 1:numel(dur)
                window = [bout.ZTstart - dur(idur), bout.ZTstart];
                windowedStateSeries = stateSeries.getWindow(time.ZT(window));
                stateRatios = windowedStateSeries.getStateRatios;
                tbl=stateRatios.State(:,{'A-WAKE','Q-WAKE','SWS','REM'});
                tbl.("A-WAKE")=seconds(tbl.("A-WAKE"));
                tbl.("Q-WAKE")=seconds(tbl.("Q-WAKE"));
                tbl.("SWS")=seconds(tbl.("SWS"));
                tbl.("REM")=seconds(tbl.("REM"));
                dd=minutes(dur(idur));
                tbl=renamevars(tbl,["A-WAKE" "Q-WAKE" "SWS" "REM"], ...
                    {sprintf('A-WAKE %dmins',dd), ...
                    sprintf('Q-WAKE %dmins',dd), ...
                    sprintf('SWS %dmins',dd), ...
                    sprintf('REM %dmins',dd), ...
                    });
                sleepFrac = stateRatios.getSleepFractionPreceeding;
                sleepFracTbl = array2table(sleepFrac.sleepFraction, ...
                    'VariableNames', sprintf("dur%dmins", minutes(dur(idur))));
                timeSeriesTbl = [timeSeriesTbl, sleepFracTbl, tbl];
            end
            boutData = [bout, timeSeriesTbl];
        end
        function boutData = processBoutAlt(~, bout, stateSeries, dur, interest,rest)
            % Processes a bout of data to extract relevant information.
            timeSeriesTbl = table();
            for idur = 1:numel(dur)
                window = [bout.ZTstart - dur(idur), bout.ZTstart];
                windowedStateSeries = stateSeries.getWindow(time.ZT(window));
                stateRatios = windowedStateSeries.getStateRatios;
                for ist=1:numel(interest)
                    st=interest(ist);
                    interestFraction = stateRatios.getStateFractionPreceeding(st,[rest interest(1:end~=ist)]);
                    interestFracTbl = array2table(interestFraction.fraction, ...
                        'VariableNames', sprintf("dur%dmins%s", ...
                        minutes(dur(idur)),st));
                    timeSeriesTbl = [timeSeriesTbl, interestFracTbl];
                end
            end
            boutData = [bout, timeSeriesTbl];
        end
    end
end
