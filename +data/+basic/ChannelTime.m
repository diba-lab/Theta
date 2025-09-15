classdef ChannelTime
    %CHANNELTIME Class to handle time series data for multiple channels
    %   Stores channel labels, time vector, and data matrix. Provides
    %   methods for plotting and processing the data.
    
    properties
        Channels    % Channel labels (categorical or other)
        Time        % Time vector (numeric or duration)
        Data        % Data matrix (channels x time)
    end
    
    methods
        function obj = ChannelTime(channel,time,data)
            %CHANNELTIME Construct an instance of this class
            %   channel: channel labels (can be empty)
            %   time: time vector
            %   data: data matrix (channels x time)
            if nargin>0
                if isempty(channel)
                    % If no channel labels provided, use default numeric labels
                    obj.Channels = categorical(1:size(data,1));
                else
                    obj.Channels = channel;
                end
                obj.Time = time;
                obj.Data = data;
            end
        end
        
        function [p] = plot(obj,color)
            %PLOT Plot the data with optional color
            %   color: RGB vector for plot color (optional)
            if nargin>1
                % If color is provided, adjust color for channels
                p1=obj.plotChannels(color+(1-color)/2);
            else
                p1=obj.plotChannels;
            end
            hold on
            % Plot error bars
            p=obj.plotErrorBar;
            if nargin>1
                % Set color for error bar plot
                p=p.setColor(color);
            end
        end

        function [p] = plotChannels(obj,color,chno)
            %PLOTCHANNELS Plot data for selected channels
            %   color: RGB vector for plot color (optional)
            %   chno: channel indices or labels (optional)
            if nargin>2
                if isnumeric(chno)
                    % Select channels by numeric index
                    obj.Data=obj.Data(chno,:);
                else
                    % Select channels by label
                    obj.Data=obj.Data(ismember(obj.Channels,chno),:);
                end
            end
            if isduration(obj.Time)
                % Convert duration to hours if needed
                obj.Time=hours(obj.Time);
            end
            if nargin>1
                % Plot with specified color
                p=plot(obj.Time,obj.Data,Color=color);
            else
                % Plot with default color
                p=plot(obj.Time,obj.Data);
            end
        end

        function [p] = plotErrorBar(obj,color)
            %PLOTERRORBAR Plot shaded error bars for the data
            %   color: RGB vector for plot color (optional)
            tbl=plotdata.ErrorBarTable(obj.Time,obj.Data);
            p=tbl.plotShaded;
            if nargin>1
                % Set color for shaded error bar
                p=p.setColor(color);
            end
        end

        function sr = getSampleRate(obj)
            %GETSAMPLERATE Calculate sample rate from time vector
            sr=1/median(diff(seconds(obj.Time)));
        end

        function obj = fillmissing(obj,n)
            %FILLMISSING Fill missing data using k-nearest neighbors
            %   n: window size (can be duration or numeric)
            if isduration(n)
                n=seconds(n)*obj.getSampleRate;
            end
            obj.Data=fillmissing(obj.Data,"knn",n);
        end

        function obj = getMedianFiltered(obj,n)
            %GETMEDIANFILTERED Apply moving median filter to data
            %   n: window size (can be duration or numeric)
            if isduration(n)
                n=seconds(n)*obj.getSampleRate;
            end
            obj.Data=smoothdata(obj.Data,2,"movmedian",n,"omitmissing");
        end

        function obj = getGaussianFiltered(obj,n)
            %GETGAUSSIANFILTERED Apply Gaussian filter to data
            %   n: window size (can be duration or numeric)
            if isduration(n)
                n=seconds(n)*obj.getSampleRate;
            end
            obj.Data=smoothdata(obj.Data,2,"gaussian",n,"omitmissing");
        end

        function obj = getHighpassFiltered(obj,freq)
            %GETHIGHPASSFILTERED Apply high-pass filter to data
            %   freq: cutoff frequency (Hz)
            obj.Data=ft_preproc_highpassfilter(obj.Data, ...
                obj.getSampleRate,freq);
        end

        function obj = getDetrended(obj)
            %GETDETRENDED Remove linear trend from data
            obj.Data=detrend(obj.Data);
        end

    end
end
