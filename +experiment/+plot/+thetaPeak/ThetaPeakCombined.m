classdef ThetaPeakCombined
    %THETAPEAKCOMBINED Combines and manages multiple ThetaPeak objects for analysis and plotting.
    %   This class holds a list of ThetaPeak objects and provides methods to
    %   merge, compare, and plot their properties in various ways.
    
    properties
        thpkList    % List of ThetaPeak objects (CellArrayList)
        Info        % Additional information (structure or object)
    end
    
    methods
        function obj = ThetaPeakCombined(thpk)
            % Constructor: Initializes the ThetaPeakCombined object.
            %   thpk: a ThetaPeak object to add to the list initially.
            obj.thpkList=CellArrayList();
            try
                obj.thpkList.add(thpk);
            catch
                % Ignore errors if thpk is not valid
            end
        end
        
        function obj = plus(obj,thpk)
            % Overloads the + operator to add a ThetaPeak object to the list.
            obj.thpkList.add(thpk);
        end

        function obj = add(obj,thpk,num)
            % Adds a ThetaPeak object to the list at a specific position.
            %   thpk: ThetaPeak object to add
            %   num: position in the list
            obj.thpkList.add(thpk,num);
        end

        function newthpks = merge(obj, thpks)
            % Merges another ThetaPeakCombined object or all objects in the list.
            %   thpks: another ThetaPeakCombined object (optional)
            if exist('thpks','var')
                newthpks=experiment.plot.thetaPeak.ThetaPeakCombined;
                if isa(thpks,'experiment.plot.thetaPeak.ThetaPeakCombined')
                    % Merge corresponding ThetaPeak objects from both lists
                    for il=1:max(obj.thpkList.length,thpks.thpkList.length)
                        thpk1=obj.thpkList.get(il);
                        thpk2=thpks.thpkList.get(il);
                        try
                            thpkmsum=thpk1.merge(thpk2);
                        catch
                            thpkmsum=thpk2.merge(thpk1);
                        end
                        newthpks=newthpks.add(thpkmsum,il);
                    end
                else
                    newthpks=obj;
                end
            else
                % Merge all ThetaPeak objects in the current list
                for ith=1:obj.thpkList.length
                    a=obj.thpkList.get(ith);
                    if ith==1
                        newthpks=a;
                    else
                        newthpks=newthpks.merge(a);
                    end
                end
            end
        end

        function cmp = compare(obj,thpks)
            % Compares ThetaPeak objects in this and another ThetaPeakCombined.
            %   thpks: another ThetaPeakCombined object
            if isa(thpks,'experiment.plot.thetaPeak.ThetaPeakCombined')
                for il=1:max(obj.thpkList.length,thpks.thpkList.length)
                    thpk1=obj.thpkList.get(il);
                    thpk2=thpks.thpkList.get(il);
                    try
                        cmp(il)=thpk1.compare(thpk2);
                    catch
                        cmp(il)=thpk2.compare(thpk1);
                    end
                end
            else
                newthpks=obj;
            end
        end

        function axsr=plotCF(obj,rows,row,col)
            % Plots the cross-frequency (CF) analysis for each ThetaPeak object.
            %   rows, row, col: subplot grid parameters or axes handles
            if ~exist('rows','var')
                rows=1;
            else
                if isa(rows,'matlab.graphics.axis.Axes')
                    axs=rows;
                end
            end
            if ~exist('row','var')
                row=1;
            end
            list=obj.thpkList;
            for isub=1:list.length
                if exist('axs','var')
                    axes(axs(isub)); %#ok<LAXES>
                    hold on
                else
                    if exist('col','var')
                        subplot(rows, col, (row-1)*col+ isub)
                    else
                        subplot(rows, list.length, (row-1)*list.length + isub)
                    end
                end
                thesub=list.get(isub);
                if ~isempty(thesub.Signal)
                    thesub.plotCF
                end
                if isub>1
                    xlabel('');
                    xticklabels([]);
                else
                    try
                        t=text(5 ,0 , obj.Info.Session.toString, ...
                            "Interpreter","none",'VerticalAlignment','bottom');
                        t.FontSize=7;
                    catch
                    end
                end
                if row<rows
                    ylabel('');
                end
                yticks([]);
                axsr(isub)=gca;
            end            
        end

        function axsr=plotSpeed(obj,rows,row,col)
            % Plots the speed for each ThetaPeak object.
            %   rows, row, col: subplot grid parameters or axes handles
            if ~exist('rows','var')
                rows=1;
            else
                if isa(rows,'matlab.graphics.axis.Axes');
                    axs=rows;
                end
            end
            if ~exist('row','var')
                row=1;
            end
            list=obj.thpkList;
            for isub=1:list.length
                if exist('axs','var')
                    axes(axs(isub)); %#ok<LAXES>
                    hold on
                else
                    if exist('col','var')
                        subplot(rows, col, (row-1)*col+ isub)
                    else
                        subplot(rows, list.length, (row-1)*list.length + isub)
                    end
                end
                thesub=list.get(isub);
                if ~isempty(thesub.Speed)
                    thesub.plotSpeed
                end
                if isub>1
                    xlabel('');
                    xticklabels([]);
                else
                    try
                        t=text(5 ,0 , obj.Info.Session.toString, ...
                            "Interpreter","none",'VerticalAlignment','bottom');
                        t.FontSize=7;
                    catch
                    end
                end
                if row<rows
                    ylabel('');
                end
                yticks([]);
                axsr(isub)=gca;
            end            
        end

        function axsr=plotCF3(obj,rows,row,col)
            % Plots a 3rd variant of cross-frequency analysis for each ThetaPeak object.
            %   rows, row, col: subplot grid parameters or axes handles
            if ~exist('rows','var')
                rows=1; 
            else
                if isa(rows,'matlab.graphics.axis.Axes');
                    axs=rows;
                end
            end
            if ~exist('row','var')
                row=1; 
            end
            list=obj.thpkList;
            for isub=1:list.length
                if exist('axs','var')
                    axes(axs(isub)); %#ok<LAXES>
                    hold on
                else
                    if exist('col','var')
                        subplot(rows, col, (row-1)*col+ isub)
                    else
                        subplot(rows, list.length, (row-1)*list.length + isub)
                    end
                end
                thesub=list.get(isub);
                if ~isempty(thesub.Signal)
                    thesub.plotCF3
                end
                if isub>1
                    xlabel('');
                    xticklabels([]);
                else
                    try
                        t=text(5 ,0 , obj.Info.Session.toString, ...
                            "Interpreter","none",'VerticalAlignment','bottom');
                        t.FontSize=7;
                    catch
                    end
                end
                if row<rows
                    ylabel('');
                end
                yticks([]);
                axsr(isub)=gca;
            end            
        end

        function axsr=plotPW(obj,rows,row,col)
            % Plots the power (PW) for each ThetaPeak object.
            %   rows, row, col: subplot grid parameters or axes handles
            if ~exist('rows','var')
                rows=1; 
            else
                if isa(rows,'matlab.graphics.axis.Axes');
                    axs=rows;
                end
            end
            
            if ~exist('row','var')
                row=1;
            end
            list=obj.thpkList;
            for isub=1:list.length
                if exist('axs','var')
                    axes(axs(isub)); %#ok<LAXES>
                    hold on
                else
                    if exist('col','var')
                        subplot(rows, col, (row-1)*col+ isub)
                    else
                        subplot(rows, list.length, (row-1)*list.length + isub)
                    end
                end
                thesub=list.get(isub);
                if ~isempty(thesub.Signal)
                    thesub.plotPW
                end
                if isub>1
                    xlabel('');
                    xticklabels([]);
                else
                    try
                        t=text(5 ,0 , obj.Info.Session.toString, ...
                            "Interpreter","none",'VerticalAlignment','bottom');
                        t.FontSize=7;
                    catch
                    end
                end
                if row<rows
                    ylabel('');
                end
                yticks([]);
                axsr(isub)=gca;
            end            
        end

        function axsr=plotPW3(obj,rows,row,col)
            % Plots a 3rd variant of power (PW) for each ThetaPeak object.
            %   rows, row, col: subplot grid parameters or axes handles
            if ~exist('rows','var')
                rows=1; 
            else
                if isa(rows,'matlab.graphics.axis.Axes');
                    axs=rows;
                end
            end
            
            if ~exist('row','var')
                row=1;
            end
            list=obj.thpkList;
            for isub=1:list.length
                if exist('axs','var')
                    axes(axs(isub)); %#ok<LAXES>
                    hold on
                else
                    if exist('col','var')
                        subplot(rows, col, (row-1)*col+ isub)
                    else
                        subplot(rows, list.length, (row-1)*list.length + isub)
                    end
                end
                thesub=list.get(isub);
                if ~isempty(thesub.Signal)
                    thesub.plotPW3
                end
                if isub>1
                    xlabel('');
                    xticklabels([]);
                else
                    try
                        t=text(5 ,0 , obj.Info.Session.toString, ...
                            "Interpreter","none",'VerticalAlignment','bottom');
                        t.FontSize=7;
                    catch
                    end
                end
                if row<rows
                    ylabel('');
                end
                yticks([]);
                axsr(isub)=gca;
            end            
        end

        function tableall=plotDurationFrequency(obj,ax,colors)
            % Plots duration vs frequency for each ThetaPeak object and returns a table.
            %   ax: axes handle for plotting
            %   colors: color array for each plot
            if ~exist('ax','var')
                ax=gca;
            end
            list=obj.thpkList;
            for isub=1:list.length
                fooof1=obj.Info.fooof(isub);
                thesub=list.get(isub);
                if ~isempty(thesub.Signal)
                    thesub.fooof.Info.episode=fooof1.Info.episode;
                    table1=thesub.plotDurationFrequency(ax,colors(isub,:));
                    if exist('tableall','var')
                        tableall=[tableall;table1];
                    else
                        tableall=table1;
                    end
                end
            end
        end

        function tableall=plotDurationFrequency3(obj,ax,colors)
            % Plots a 3rd variant of duration vs frequency for each ThetaPeak object and returns a table.
            %   ax: axes handle for plotting
            %   colors: color array for each plot
            if ~exist('ax','var')
                ax=gca;
            end
            list=obj.thpkList;
            for isub=1:list.length
                fooof1=obj.Info.fooof(isub);
                thesub=list.get(isub);
                if ~isempty(thesub.Signal)
                    thesub.fooof.Info.episode=fooof1.Info.episode;
                    table1=thesub.plotDurationFrequency3(ax,colors(isub,:));
                    if exist('tableall','var')
                        tableall=[tableall;table1];
                    else
                        tableall=table1;
                    end
                end
            end
        end

    end
end
