classdef MSPlot_ < matlab.mixin.Copyable
    %MSPLOT_ Provides graphical capabilities to visualize data from a MSModel_
    %   Plot Types:
    %           - MH plots (specify MH pairs, all plotted on same axis)
    %
    
    properties
        fig_num         %number of figure to plot in
        plot_handle     %can be returned from plot functions
        frame_data      %used to create animations
        animation_view  %view for animations (2,3 or [Az El]), see: view
        file_name       %name to output file (image or movie)
        save_frames     %save frames to output a movie file (true / false)
    end
    properties (SetAccess = protected)
        Source_field_names
        State_names
    end
    
    properties (Access = protected)
        model              
        animation_zoom  %zoom level for animations
    end
    
    methods
        %% Default constructor
        function obj = MSPlot_(model)
            %This class just stores various plotting functions
            obj.model = model; %store handle to model
            obj.fig_num = 1; %default figure number
            
            %get source field names
            sources = properties( obj.model.particles{1}.SourceFields);
            count = 0;
            for i = 1:numel(sources)
                temp = obj.model.particles{1}.SourceFields.(sources{i});
                for j = 1:numel(temp.Property_List)
                    count = count+1;
                    list{count} = [sources{i},'.',temp.Property_List{j}];                    
                end
            end
            obj.Source_field_names = list;
            
            %get names of State fields (m1, m2, etc.)
            states = properties(obj.model.particles{1}.State);
            states = setdiff(states, {'m1','m2','m3','Property_List'});
            obj.State_names ={states{:}};
            
            %default options
            obj.save_frames = false;
            obj.animation_view = 3;
            obj.animation_zoom = 1; 
            obj.file_name = 'file';
        end
        
        %% Plot initial geometry / prepare animation plot
        function [obj,h] = Plot_Initial_State(obj,particle_nums,varargin)
         
            %initial input check
            switch isnumeric(particle_nums)
                case 0
                    str = sprintf('The first input must be a 1 dimensional array of numbers indicating which particles to plot.');
                    warning(str)
            end             
            [Prop,Val] = MS.PropertyValue(varargin); %get extra input arguments
            %check input properties
            
            %Load relevant data
            m1 = cell(1,numel(particle_nums));
            m2 = m1; m3 = m1; X = m1; dims = m1; shape = m1;
            for i = 1:numel(particle_nums)
                particle = obj.model.particles{particle_nums(i)};
                m1{i} = particle.State.m1;
                m2{i} = particle.State.m2;
                m3{i} = particle.State.m3;
                X{i} = particle.Properties.Location;
                dims{i} = particle.Properties.Dims;
                shape{i} = particle.Properties.Shape;                
            end
            
            
            %% Initialize plot window
            %setup figure
            h.fig = figure(obj.fig_num);
            clf
            h.ax = gca;
            hold on
            %plot initial shapes / lines / points
            unit_scale = cell(1,numel(m1)); %scale for m vector
            for i = 1:numel(m1)                
                %initial shapes
                switch shape{i}
                    case 'Rectangle'                        
                        [x,y,z] = myRectangle('Position',[x y z w t h]);
                    case 'Ellipse'
                        [x,y,z] = ellipsoid(...
                            X{i}(1),X{i}(2),X{i}(3),...
                            dims{i}(1)/2, dims{i}(2)/2, dims{i}(3)/2, 40);
                end
                h.shapes{i} = surf(x,y,z);
                h.shapes{i}.EdgeAlpha = 0.05;
                %                 h.shapes{i}.FaceColor = [0.3 0.75 0.9];
                %                 h.shapes{i}.FaceAlpha = 0.2;
                h.shapes{i}.FaceColor = 'w';
                h.shapes{i}.FaceAlpha = 0.01;
                hold on
                
                %plot unit sphere around shape
                unit_scale{i} = max([dims{:}])/2;
                [x,y,z] = ellipsoid(...
                    X{i}(1),X{i}(2),X{i}(3),...
                    unit_scale{i}, unit_scale{i}, unit_scale{i});
                h.unit_sphere{i} = surf(x,y,z);
                h.unit_sphere{i}.EdgeAlpha = 0.15;
                h.unit_sphere{i}.FaceColor = 'w';
                h.unit_sphere{i}.FaceAlpha = 0.03;
                
                %initial magnetization
                [h.arrowHead{i},h.arrowStem{i}] =...
                    quiver3D([X{i}(1),X{i}(2),X{i}(3)],... %from mathworks file exchange
                    [m1{i}(1),m2{i}(1),m3{i}(1)]*unit_scale{i},...
                    'r');
                
                %trace for each arrow                
                h.trace_line{i} = line(X{i}(1),X{i}(2),X{i}(3),'LineStyle','-','Color','b','LineWidth',1.25);                  %ALLOW USER TO EDIT
                h.trace_tip{i} = line(X{i}(1),X{i}(2),X{i}(3),'Color','r','Marker','*','MarkerSize',7);       %ALLOW USER TO EDIT
                                
            end
            %text to update 
            h.text = text(0.15,0.25,0.9,'Initial State ', 'Units','Normalized');
            h.unit_scale =unit_scale; % save for m vector / unit sphere;
            axis square
            lighting phong
            camlight head
            
            
            %make each axis scale the same
%             axis equal
            [max_scale,max_ind] = max([diff(h.ax.XLim), diff(h.ax.YLim), diff(h.ax.ZLim)]);
            max_scale = max_scale/2;
            scale_factor = 0.9;                                                             %ALLOW USER TO EDIT
            switch max_ind
                case 1 %x axis is currently largest
                    h.ax.YLim = mean(h.ax.YLim) + scale_factor*[-max_scale, max_scale];
                    h.ax.ZLim = mean(h.ax.ZLim) + scale_factor*[-max_scale, max_scale];
                case 2 %y axis is currently largest
                    h.ax.XLim = mean(h.ax.XLim) + scale_factor*[-max_scale, max_scale];
                    h.ax.ZLim = mean(h.ax.ZLim) + scale_factor*[-max_scale, max_scale];
                case 3 %z axis is currently largest
                    h.ax.XLim = mean(h.ax.XLim) + scale_factor*[-max_scale, max_scale];
                    h.ax.YLim = mean(h.ax.YLim) + scale_factor*[-max_scale, max_scale];
            end
            
            %get rid of 3D box
            h.ax.Visible = 'off';
            
            %add custom planes
%             height = mean(h.ax.ZLim);
%             h.ground_plane = patch(...
%                 [h.ax.XLim(1) h.ax.XLim(2) h.ax.XLim(2) h.ax.XLim(1)],...
%                 [h.ax.YLim(1) h.ax.YLim(1) h.ax.YLim(2) h.ax.YLim(2)],...
%                 [height height height height],'k');
%             h.ground_plane.FaceColor = [.83 .82 .78];
%             h.ground_plane.FaceAlpha = .4;
            
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            
            
            grid on
            obj.plot_handle = h;
            
            zoom(obj.animation_zoom)
            view(obj.animation_view)
            
            drawnow
            
        end %Plot_Initial_State
        
        %% Plot data for single particles
        function obj = Particle_Plot(obj,varargin)
            %M vs H, M vs time, Angle vs H, Angele vs time
            % data must be entered as a list of triplets
            % 1 - particle number
            % 2 - X data
            % 3 - Y data
            
            %initial input check
            switch mod(numel(varargin),3)
                case 0
                    %get input triplets
                    [particle_num, Xstr, Ystr] = MS.PropertyValue3(varargin);
                otherwise
                    msg = sprintf(['Plot_Particle accepts lists containing sequenes of three items:',...
                        '\n 1 - Particle Number \n 2 - X Data (string) \n 3 - Y Data (string)',...
                        '\nMultiple triplets may be entered, triplets with the same X data will be plotted together']);
                    error(msg)
                    
            end
            
            % Sort inputs based on X data (so all plots with the same x
            % data are plotted together
            [Xstr,sort_ind] = sort(Xstr);
            Ystr = Ystr(sort_ind);
            particle_num = particle_num(sort_ind);
            
            %% Plot results for each triplet
            fig = gobjects(1,numel(unique(Xstr)));%initialize figure handles
            ax = fig;   %initialize axis handles
            leg = ax;   %initialize legend handles
            h = gobjects(1,numel(particle_num)); %initialize line handles
            
            %temp variables
            count1 = 1; %figure number
            
            %create first figure object
            fig(count1) = figure(obj.fig_num);
            clf
            ax(count1) = gca;
            xlabel(Xstr{1})
            
            for i = 1:numel(particle_num)
                %load particle
                particle = obj.model.particles{particle_num{i}}; %#ok<NASGU>
                %see if the X data is a match from the previous plot
                switch ( i > 1 ) && ( ~strcmp(Xstr{i} , Xstr{i-1}) )
                    case 1 %if not, create a new figure
                        %first show the legend / grid
                        legend('show')
                        leg(count1) = findobj(fig(count1).Children,'Tag','legend');
                        grid on
                        %make new figure
                        count1 = count1 + 1;
                        fig(count1) = figure(obj.fig_num+count1-1);
                        clf
                        ax(count1) = gca;
                        xlabel(Xstr{i})
                    case 0
                        if i > 1
                            ax(count1).ColorOrderIndex = mod(ax(count1).ColorOrderIndex,size(ax(count1).ColorOrder,1)) +1; %advance color
                        end
                        hold all;
                end
                %get X values
                if ismember(Xstr{i},obj.Source_field_names)
                    X{i} = eval(['particle.SourceFields.',Xstr{i}]);
                elseif ismember(Xstr{i},obj.State_names)
                    X{i} = eval(['particle.State.',Xstr{i}]);
                end
                %get Y values
                if ismember(Ystr{i},obj.Source_field_names)
                    Y{i} = eval(['particle.SourceFields.',Ystr{i}]);
                elseif ismember(Ystr{i},obj.State_names)
                    Y{i} = eval(['particle.State.',Ystr{i}]);
                end
                
                h(i) = line(X{i},Y{i});
                h(i).Color = ax(count1).ColorOrder(ax(count1).ColorOrderIndex,:);
                h(i).DisplayName = strrep(['P',num2str(particle_num{i}),': ',Ystr{i}],'_',' ');
                
                %show legend for final plot
                if i == numel(particle_num)
                    legend('show')
                    leg(count1) = findobj(fig(count1).Children,'Tag','legend');
                    grid on
                end
            end
            
            %return handles
            obj.plot_handle.fig = fig;
            obj.plot_handle.ax = ax;
            obj.plot_handle.lines = h;
            obj.plot_handle.legend = leg;
            
        end %Particle_Plot()
        
        %% Animate all spins
        function [obj,h] = Animate_Spins(obj,particle_nums,data_type,varargin)
            %Plot m vs time for the indicated particles
            % data should be entered as an array containing the numbers of
            % which particles to plot, as well as a string indicating if
            % the static or dynamic values should be plotted            
            
            %% initial setup
            %initial input check
            switch isnumeric(particle_nums) && all(ismember(data_type,{'static','dynamic'})) && (numel(data_type)==numel(particle_nums) || numel(data_type)==1 )
                case 0
                    str = sprintf(['The first input must be a 1 dimensional array of numbers indicating which particles to plot.',...
                        '\nThe second input must be a cell array of strings with ''static'' or ''dynamic''. The cell array can have\n',...
                        'either one element, or as many elements as the number array, so you can combine static / dynamic plots)']);
                    warning(str)
            end
            %make sure data_type is a cell with the same number of entries
            %as particle_nums. If the initial size is 1, duplicate the same
            %entry for all particles
            if ~iscell(data_type); data_type = {data_type}; end
            if numel(data_type)==1; data_type = repmat(data_type,size(particle_nums)); end
            
            [Prop,Val] = MS.PropertyValue(varargin); %get extra input arguments
            %check input properties
            
            %Load relevant data
            m1 = cell(1,numel(particle_nums));
            m2 = m1; m3 = m1; X = m1; dims = m1; shape = m1;
            for i = 1:numel(particle_nums)
                particle = obj.model.particles{particle_nums(i)};
                switch data_type{i}
                    case 'static'
                        m1{i} = particle.State.Stat_m1;
                        m2{i} = particle.State.Stat_m2;
                        m3{i} = particle.State.Stat_m3;
                    case 'dynamic'
                        m1{i} = particle.State.Dyn_m1;
                        m2{i} = particle.State.Dyn_m2;
                        m3{i} = particle.State.Dyn_m3;
                end
                X{i} = particle.Properties.Location;
                dims{i} = particle.Properties.Dims;
                shape{i} = particle.Properties.Shape;
                switch i > 1
                    case 1
                        size_check = numel(m1{i}) ~= numel(m1{i-1});
                        if size_check;
                            error('When mixing static and dynamic plots, the source fields must have the same number of elements');
                        end
                end
            end
            
            
            %% Initialize plot window
            [~,h] = Plot_Initial_State(obj,particle_nums,data_type,varargin);
            unit_scale = h.unit_scale;
            
            %% Animate            
            npts = numel(m1{1});
            nparts = numel(m1);
                        
            switch npts >= 1000 
                case 1
                    spacing = round(npts / 1000); %plot a total of 1000 points
%                     spacing = 1;
                case 0
                    spacing = 1;
            end
            
            time_points = 1:spacing:npts;
            xtip = cell(1,nparts);
            count = 1;
            for i = 1:nparts; xtip{i} = zeros(numel(time_points),3); end;
            
            %initialize frame structure
            if obj.save_frames
                obj.frame_data = getframe(gca);    
                obj.frame_data(numel(time_points)) = getframe(gca);
            end
                
            for i = time_points         %loop over all time points
                for j = 1:nparts    %loop over particles
                    %update arrow
                    [ Xstem,Xhead ] = arrow3Dupdate(...
                        [X{j}(1),X{j}(2),X{j}(3)],...
                        [m1{j}(i),m2{j}(i),m3{j}(i)]*unit_scale{j} );
                    
                    h.arrowHead{j}.XData = Xhead{1}{1};
                    h.arrowHead{j}.YData = Xhead{1}{2};
                    h.arrowHead{j}.ZData = Xhead{1}{3};
                    
                    h.arrowStem{j}.XData = Xstem{1}{1};
                    h.arrowStem{j}.YData = Xstem{1}{2};
                    h.arrowStem{j}.ZData = Xstem{1}{3};
                    
                    %update trace
                    xtip{j}(count,:) = [X{j}(1),X{j}(2),X{j}(3)] + ...
                        [m1{j}(i),m2{j}(i),m3{j}(i)]*unit_scale{j};
                    h.trace_line{j}.XData = xtip{j}(1:count,1);
                    h.trace_line{j}.YData = xtip{j}(1:count,2);
                    h.trace_line{j}.ZData = xtip{j}(1:count,3);
                    
                    h.trace_tip{j}.XData = xtip{j}(count,1);
                    h.trace_tip{j}.YData = xtip{j}(count,2);
                    h.trace_tip{j}.ZData = xtip{j}(count,3);                                                            
                end
                %update text
                percent_done = round(i/max(time_points)*100);
                h.text.String = sprintf('Completed: %d%%',percent_done);
                
                %Save frames to create movie
                if obj.save_frames
                    obj.frame_data(count) = getframe(gca);
                end
                drawnow
                count = count + 1;
            end
            
        end %Animate_Spins
                  
        %% Energy plots
           function [obj,h] = animate_energy_2d(obj,particle_num,data_type,varargin)
            %Plot the 2D energy surface for the indicated particles
            % data should be entered as an array containing the numbers of
            % which particles to plot
            
            %% input check
            %initial input check
            switch isnumeric(particle_num) && all(ismember(data_type,{'static','dynamic'})) && (numel(data_type)==numel(particle_num) || numel(data_type)==1 )
                case 0
                    str = sprintf(['The first input must be a 1 dimensional array of numbers indicating which particles to plot.',...
                        '\nThe second input must be a cell array of strings with ''static'' or ''dynamic''. The cell array can have\n',...
                        'either one element, or as many elements as the number array, so you can combine static / dynamic plots)']);
                    error(str)
            end
            if numel(particle_num)>1
                error('Only one particle at a time can be analyzed using animate_energy_2d')
            end
            
            %make sure data_type is a cell with the same number of entries
            %as particle_nums. If the initial size is 1, duplicate the same
            %entry for all particles
            if ~iscell(data_type); data_type = {data_type}; end
            if numel(data_type)==1; data_type = repmat(data_type,size(particle_num)); end
            
            [Prop,Val] = MS.PropertyValue(varargin); %get extra input arguments
            %check input properties
            
            %% Load relevant data           
            particle = obj.model.particles{particle_num};
            switch data_type{1}
                case 'static'
                    H1 = particle.SourceFields.static.H1;
                    H2 = particle.SourceFields.static.H2;
                    H3 = particle.SourceFields.static.H3;
                    S1 = particle.SourceFields.static.S1;
                    S2 = particle.SourceFields.static.S2;
                    S3 = particle.SourceFields.static.S3;
                    S4 = particle.SourceFields.static.S4;
                    S5 = particle.SourceFields.static.S5;
                    S6 = particle.SourceFields.static.S6;
                    Sp1 = particle.SourceFields.static.Sp1;
                    Sp2 = particle.SourceFields.static.Sp2;
                    Sp3 = particle.SourceFields.static.Sp3;
                    p1 = particle.SourceFields.static.p1;
                    p2 = particle.SourceFields.static.p2;
                    p3 = particle.SourceFields.static.p3;
                case 'dynamic'
                    H1 = particle.SourceFields.dynamic.H1;
                    H2 = particle.SourceFields.dynamic.H2;
                    H3 = particle.SourceFields.dynamic.H3;
                    S1 = particle.SourceFields.dynamic.S1;
                    S2 = particle.SourceFields.dynamic.S2;
                    S3 = particle.SourceFields.dynamic.S3;
                    S4 = particle.SourceFields.dynamic.S4;
                    S5 = particle.SourceFields.dynamic.S5;
                    S6 = particle.SourceFields.dynamic.S6;
                    Sp1 = particle.SourceFields.dynamic.Sp1;
                    Sp2 = particle.SourceFields.dynamic.Sp2;
                    Sp3 = particle.SourceFields.dynamic.Sp3;
                    p1 = particle.SourceFields.dynamic.p1;
                    p2 = particle.SourceFields.dynamic.p2;
                    p3 = particle.SourceFields.dynamic.p3;
            end
            
            %% grid out m vector to plot
            npts = 1e3;
            theta = linspace(0,360,npts);
            m1 = cosd(theta);
            m2 = sind(theta);
            m3 = zeros(size(m1));
            
            %% Compute Energy Terms
            U_H = @(i) particle.Energy.U_Zeeman(m1,m2,m3,H1(i),H2(i),H3(i));
            U_D = @(i) particle.Energy.U_Demag(m1,m2,m3);
            U_ME = @(i) particle.Energy.U_ME(m1,m2,m3,...
               S1(i),S2(i),S3(i),S4(i),S5(i),S6(i));
           switch nargin(particle.Energy.U_EB) 
               case 3
                   U_EB = @(i) particle.Energy.U_EB(m1,m2,m3);                   
               case 5
                   U_EB = @(i) particle.Energy.U_EB(m1,m2,m3,S1(i),S2(i));
               case 7
                   U_EB = @(i) particle.Energy.U_EB(m1,m2,m3,Sp1(i),Sp3(i),p1(:,i),p3(:,i));
               case 9
                   U_EB = @(i) particle.Energy.U_EB(m1,m2,m3,S1(i),S2(i),S3(i),S4(i),S5(i),S6(i));
           end
            U_UNI = @(i) particle.Energy.U_uni(m1,m2,m3);       
            U_D_UNI = @(i) particle.Energy.U_Demag(m1,m2,m3) + particle.Energy.U_uni(m1,m2,m3);       
            U_MCA = @(i) particle.Energy.U_MCA(m1,m2,m3); 
           switch nargin(particle.Energy.U_EB) 
               case 7
                   U_TOT = @(i) particle.Energy.U_total(H1(i),H2(i),H3(i),m1,m2,m3,...
                       p1(:,i),p3(:,i), S1(i),S2(i),S3(i),S4(i),S5(i),S6(i),Sp1(i),Sp3(i));
               otherwise
                   U_TOT = @(i) particle.Energy.U_total(H1(i),H2(i),H3(i),m1,m2,m3,...
                       S1(i),S2(i),S3(i),S4(i),S5(i),S6(i));
           end
            
            %% Determine energy bounds
%             U_min = zeros(1,numel(H1)); U_max = U_min; 
            for i = 1:numel(H1)
                U_min(i,:) = min([U_H(i)',U_D(i)',U_ME(i)',U_EB(i)',U_UNI(i)',U_D_UNI(i)',U_MCA(i)',U_TOT(i)']);                
                U_max(i,:) = max([U_H(i)',U_D(i)',U_ME(i)',U_EB(i)',U_UNI(i)',U_D_UNI(i)',U_MCA(i)',U_TOT(i)']);
            end
            U_min = min(U_min);
            U_max = max(U_max);
            U_range = U_max - U_min;
            
            %% Initialize plot window
            fig = gobjects(1);  %initialize figure handles
            ax = fig;           %initialize axis handles
            leg = gobjects(2); 	%initialize legend handles
            
            %create first figure object
            fig = figure(obj.fig_num);
            clf
            
            ax(1) = gca;%subplot(2,1,1);
            hline{1} = plot(theta,[...
                U_ME(1)'-U_min(3),...
                U_MCA(1)'-U_min(7),...
                U_TOT(1)'-U_min(8)],'LineWidth',2);
            %             hline{1} = plot(theta,[U_H(1)'-U_min(1),...
            %                 U_D(1)'-U_min(2),...
            %                 U_ME(1)'-U_min(3),...
            %                 U_EB(1)'-U_min(4),...
            %                 U_UNI(1)'-U_min(5),...
            %                 U_D_UNI(1)'-U_min(6),...
            %                 U_MCA(1)'-U_min(7),...
            %                 U_TOT(1)'-U_min(8)]);
            %             hline{1}(1).DisplayName = 'H';
            %             hline{1}(2).DisplayName = 'D';
            %             hline{1}(3).DisplayName = 'ME';
            %             hline{1}(3).Marker = '.';
            %             hline{1}(4).DisplayName = 'EB';
            %             hline{1}(5).DisplayName = 'UNI';
            %             hline{1}(6).DisplayName = 'D+UNI';
            %             hline{1}(7).DisplayName = 'MCA';
            %             hline{1}(8).DisplayName = 'TOT';
            
            hline{1}(1).DisplayName = 'U_{STRAIN}';
            hline{1}(2).DisplayName = 'U_{CRYSTAL}';
            hline{1}(3).DisplayName = 'U_{TOTAL}';
            leg(1) = legend('show');
            
            uistack(hline{1}(1),'bottom')
            
            
            xlim([0, 360])
            ylim([0 max(U_range)])
%             ylim([10^(log10(U_range)-1) U_range])
            xlabel('\theta (deg)')
            ylabel('U (J/m^3)')
            ax.FontSize = 14;
            ax.FontWeight = 'bold';
            
            grid on
%             ax(1).YScale = 'log';
            
%             ax(2) = subplot(2,1,2);            
%             hline{2}(1) = energy_plot(U_H,theta,min(U_min),U_range(end),1,'H');
%             hold all
%             hline{2}(2) = energy_plot(U_D,theta,min(U_min),U_range(end),1,'D');
%             hline{2}(3) = energy_plot(U_ME,theta,min(U_min),U_range(end),1,'ME');
%             hline{2}(4) = energy_plot(U_EB,theta,min(U_min),U_range(end),1,'EB');
%             hline{2}(5) = energy_plot(U_UNI,theta,min(U_min),U_range(end),1,'UNI');            
%             hline{2}(6) = energy_plot(U_TOT,theta,min(U_min),U_range(end),1,'TOT');                                  
%            
%             xlim([-1,1])
%             ylim([-1 1])            
%             xlabel('m_x')
%             ylabel('m_y')
            
            %% Animate     
            npts = numel(H1);            
            
            switch npts >= 1000
                case 1
                    spacing = round(npts / 1000); %plot a total of 1000 points
                    %                     spacing = 1;
                case 0
                    spacing = 1;
            end
            
            %initialize frame structure
            if obj.save_frames
                obj.frame_data = getframe(gcf);
                obj.frame_data(floor(numel(H1)/spacing)) = getframe(gcf);
            end
            %animate frames
            count = 1;
            for j = 1:1
                for i = 1:spacing:numel(H1)
                    hline{1}(1).YData = U_ME(i)'-U_min(3);                    
                    hline{1}(2).YData = U_MCA(i)'-U_min(7);
                    hline{1}(3).YData = U_TOT(i)'-U_min(8);
                    
%                     subplot(2,1,1)
%                     hline{1}(1).YData = U_H(i)'-U_min(1);
%                     hline{1}(2).YData = U_D(i)'-U_min(2);
%                     hline{1}(3).YData = U_ME(i)'-U_min(3);
%                     hline{1}(4).YData = U_EB(i)'-U_min(4);
%                     hline{1}(5).YData = U_UNI(i)'-U_min(5);                    
%                     hline{1}(7).YData = U_MCA(i)'-U_min(7);
%                     hline{1}(8).YData = U_TOT(i)'-U_min(8);                    
%                     subplot(2,1,2)
%                     [X,Y] = pol2cart(pi/180*theta',(U_H(i)'- min(U_min))/U_range(end) );
%                     hline{2}(1).XData = X;  hline{2}(1).YData = Y;
%                     
%                     [X,Y] = pol2cart(pi/180*theta',(U_D(i)'- min(U_min))/U_range(end) );
%                     hline{2}(2).XData = X;  hline{2}(2).YData = Y;
%                     
%                     [X,Y] = pol2cart(pi/180*theta',(U_ME(i)'- min(U_min))/U_range(end) );
%                     hline{2}(3).XData = X;  hline{2}(3).YData = Y;
%                     
%                     [X,Y] = pol2cart(pi/180*theta',(U_EB(i)'- min(U_min))/U_range(end) );
%                     hline{2}(4).XData = X;  hline{2}(4).YData = Y;
%                     
%                     [X,Y] = pol2cart(pi/180*theta',(U_UNI(i)'- min(U_min))/U_range(end) );
%                     hline{2}(5).XData = X;  hline{2}(5).YData = Y;
%                     
%                     [X,Y] = pol2cart(pi/180*theta',(U_TOT(i)'- min(U_min))/U_range(end) );
%                     hline{2}(6).XData = X;  hline{2}(6).YData = Y;

                    %Save frames to create movie
                    if obj.save_frames
                        obj.frame_data(count) = getframe(gcf);
                        count = count + 1;
                    end

                    drawnow
                end
            end
            
        end
                       
        %% Save movie
        function [obj,h] = make_movie(obj,varargin)
%             obj.save_frames = true;
            obj.animation_zoom = 1;
%             obj.animation_view = [-21 22]; 
%             obj.animation_view = [3]; 
%             fname = 'bennet-clocking-2';
            %generate data
%             [~,h] = obj.Animate_Spins(particle_nums,data_type);
            
            %write movie
            frames = obj.frame_data;
            if isempty(frames)
                error('must save frames before writing a movie')
            end
            vid_profile = 'Motion JPEG AVI'; %SEE: VideoWriter - profiles 'MPEG-4', 'Motion JPEG AVI', 'Motion JPEG 2000', 'Uncompressed AVI'
            writerObj = VideoWriter(obj.file_name,vid_profile);
            writerObj.Quality = 100;
            open(writerObj)
            for i = 1:numel(frames)
                writeVideo(writerObj,frames(i));
            end
            close(writerObj)
            
            %restore default to NOT storing movie frames
            obj.save_frames = false;
        end
    
      %%
        
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               Local Functions                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h] = energy_plot(U,theta,umin,urange,i,str)
[X,Y] = pol2cart(pi/180*theta',(U(i)'-umin)/urange );
h = plot(X,Y,'DisplayName',str);

end

%% Plot Data
% %create a spherical grid of possble m orientations
% nmesh = 100;
% [THETA, PHI] = meshgrid(linspace(0,pi,nmesh/2),linspace(-pi,pi,100));
% % R = ones(size(THETA));
% % [M1,M2,M3] = sph2cart(PHI,pi/2-THETA,R);
%
% %setup debugging plot
% figure(1)
% clf
% subplot(1,2,1)
% h1 = surf(THETA*180/pi,PHI*180/pi,ones(size(THETA))); hold all
% h1.EdgeAlpha = 0.1; h1.FaceColor = 'interp';
% h2 = plot3(0,0,0,'r*');
% xlabel('\Theta')
% ylabel('\Phi')
% zlabel('Energy (J/m^3')
% xlim([0 180])
% ylim([-180 180])
% view(2)
% grid on
%
% subplot(1,2,2)
% h3 = plot3(0,0,0,'bo'); hold all
% h4 = plot3(0,0,0,'r:');
% xlabel('m1')
% ylabel('m2')
% zlabel('m3')
% xlim([-1 1])
% ylim([-1 1])
% zlim([-1 1])
% grid on
%
% for j = 1:nparts
%     for i = 1:npts
%         Uplot = @(theta,phi) particles{j}.Energy.U_total(...
%             particles{j}.SourceFields.dynamic.H1(i),...
%             particles{j}.SourceFields.dynamic.H2(i),...
%             particles{j}.SourceFields.dynamic.H3(i),...
%             M1(theta,phi),M2(theta,phi),M3(theta,phi),...
%             particles{j}.SourceFields.dynamic.S1(i),...
%             particles{j}.SourceFields.dynamic.S2(i),...
%             particles{j}.SourceFields.dynamic.S3(i),...
%             particles{j}.SourceFields.dynamic.S4(i),...
%             particles{j}.SourceFields.dynamic.S5(i),...
%             particles{j}.SourceFields.dynamic.S6(i) );
%         h1.ZData = (Uplot(THETA,PHI)-min(min(Uplot(THETA,PHI))))/...
%             max(max(abs(Uplot(THETA,PHI)-min(min(Uplot(THETA,PHI))) )));
% %         h2.XData = theta(i,j)*180/pi;
% %         h2.YData = phi(i,j)*180/pi;
% %         h2.ZData = (Uplot(theta(i,j),phi(i,j))-min(min(Uplot(theta(i,j),phi(i,j)))))/...
% %             (max(max(abs(Uplot(theta(i,j),phi(i,j))-min(min(Uplot(theta(i,j),phi(i,j)))) )))+eps);
%         h3.XData = m1(i,j);
%         h3.YData = m2(i,j);
%         h3.ZData = m3(i,j);
%         h4.XData = m1(1:i,j);
%         h4.YData = m2(1:i,j);
%         h4.ZData = m3(1:i,j);
%
%         drawnow
%     end
% end
