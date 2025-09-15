classdef TwoTableComparison
    %TWOTABLECOMPARISON Summary of this class goes here
    %   Detailed explanation goes here

    properties
        table1
        table2
    end

    methods
        function obj = TwoTableComparison(tbl1,tbl2)
            %TWOTABLECOMPARISON Construct an instance of this class
            %   Detailed explanation goes here
            obj.table1 = tbl1;
            obj.table2 = tbl2;
        end

        function [] = get(obj,num_permutations, alpha)
            condition1_data=obj.table1;
            condition2_data=obj.table2;
            time_points_of_interest=1:size(condition2_data,2);
            % Perform a cluster-based permutation test to compare two sets of ERP channels

            % Determine the number of subjects in each condition
            num_subjects_condition1 = size(condition1_data, 1);
            num_subjects_condition2 = size(condition2_data, 1);

            % Determine the minimum number of subjects in either condition
            min_num_subjects = min(num_subjects_condition1, num_subjects_condition2);

            % Initialize a variable to store cluster-level statistics
            cluster_stats = zeros(1, num_permutations);

            % Calculate the observed test statistic (e.g., t-statistic) for each time point
            observed_statistic =mean(condition1_data, 1,"omitmissing") - mean(condition2_data, 1,"omitmissing");

            % Perform the permutation test
            for perm = 1:num_permutations
                % Randomly select subjects for each condition with replacement
                condition1_indices = randperm(num_subjects_condition1, min_num_subjects);
                condition2_indices = randperm(num_subjects_condition2, min_num_subjects);

                % Extract data for the selected subjects
                condition1_sampled = condition1_data(condition1_indices, :);
                condition2_sampled = condition2_data(condition2_indices, :);

                % Calculate the test statistic for the permuted data
                permuted_statistic = mean(condition1_sampled, 1,"omitmissing") - mean(condition2_sampled, 1,"omitmissing");

                % Compute the cluster-level statistic using a threshold (e.g., t-threshold)
                threshold = 15; % Adjust this threshold based on your data and significance level
                clusters = bwconncomp(permuted_statistic > threshold);

                % Calculate cluster-level statistics (sum of test statistics within clusters)
                cluster_sum = cellfun(@(x) sum(permuted_statistic(x)), clusters.PixelIdxList);

                % Store the maximum cluster-level statistic
                if ~isempty(cluster_sum)
                    cluster_stats(perm) = max(cluster_sum);
                end
            end

            % Calculate the critical threshold for significance
            threshold_critical = prctile(cluster_stats, (1 - alpha) * 100);

            % Find clusters in the observed data that exceed the critical threshold
            observed_clusters = bwconncomp(abs(observed_statistic) > threshold);

            % Calculate cluster-level statistics for observed clusters
            observed_cluster_sum = cellfun(@(x) sum(observed_statistic(x)), observed_clusters.PixelIdxList);

            % Find significant clusters
            significant_clusters = find(abs(observed_cluster_sum) > threshold_critical);

            % Get significant time points
            significant_time_points = cellfun(@(x) x, observed_clusters.PixelIdxList(significant_clusters), 'UniformOutput', false);

            % Plot both data and test statistics overlayed
            figure;
            plot(time_points_of_interest, observed_statistic, 'k', 'LineWidth', 2);
            hold on;
            for i = 1:length(significant_clusters)
                cluster_range = significant_time_points{i};

                % Overlay test statistics (significant clusters) as a red line
                plot(time_points_of_interest(cluster_range), observed_statistic(cluster_range), 'r', 'LineWidth', 2);

                % Add asterisks to significant data points
                plot(time_points_of_interest(cluster_range), observed_statistic(cluster_range), 'r*', 'MarkerSize', 8);
            end
            xlabel('Time Points');
            ylabel('Test Statistic');
            title('Cluster-Based Permutation Test');
            legend('Observed Data', 'Significant Clusters');

            % Display the critical threshold
            fprintf('Critical Threshold: %.2f\n', threshold_critical);
        end
    end
end
