classdef MS < handle
    %MS Is a utility class that holds useful functions to be accessed
    %by multiple other classes    
    
    %     properties
    %
    %     end
    
    %     methods
    %         function obj = MSFunc_(obj) %default constructor
    %
    %         end
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                          FUNCTIONS                                %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public, Static = true)
        
        %% PropertyValue pairs / triples
        function [prop, val] = PropertyValue(input)
            %PROPERTYVALUE separate the input cell array into property-value pairs
            %   Input showd be an nx2 or 2xn cell array of property-value pairs
            
            %make sure input has even number of elements
            n = length(input);
            if mod(n,2) ~= 0
                error('Inputs must be in form of ("Property","Value") pairs')
            else
                prop = input(1:2:n-1);   %Properties
                val = input(2:2:n);      %Values
            end
            
        end
        function [particle_num, Xstr, Ystr] = PropertyValue3(input)
            %PROPERTYVALUE3 separate the input cell array into property-value triplets
            %   Input should be an 1x3^n cell array of property-value pairs
            
            %make sure input has even number of elements
            n = length(input);
            if mod(n,3) ~= 0
                error('Inputs must be in form of ("Particle","X-data","Y-data") triplets')
            else
                particle_num= input(1:3:n-2);   %Properties
                Xstr = input(2:3:n-1);      %Values
                Ystr = input(3:3:n);      %Values
            end
            
        end

        %% Create_Fields
        function [ varargout ] = Create_Field( amplitude, ave, shape, pts_per_cycle, npts )
            %CREATE_FIELD has a set of built in waveforms that can be used to create
            %the various input fields.
            %   Inputs:
            %       ave - average value of waveform
            %       amplitude - amplitude of waveformf (half peak to peak)
            %       shape - sin, cos, ramp, or const
            %       pts_per_cycle - number of data points per cycle
            %       npts - overall length of field ( note: the number of cycles is
            %       npts/pts_per_cycle)
            %   Output:
            %       n arrays of size (1 x npts). n is the length of the input amplitdue
            %       array
            %   Notes: If ave or shape don't have the same number of elements as
            %   amplitdue, then they must have only one value, which will be used for
            %   all arrays.
            
            if nargout ~= numel(amplitude)
                error('Number of outputs must be the same as the number of amplitudes');
            end
            
            %if only one shape is entered, apply for all fields
            switch numel(amplitude) ~= numel(shape) || numel(amplitude) ~= numel(ave)
                case 1
                    switch numel(shape)
                        case 1
                            shape = repmat(shape,numel(amplitude),1);
                        otherwise
                            error('Number of input shapes must equal number of input amplitudes, or be equal to 1')
                    end
                    
                    switch numel(ave)
                        case 1
                            ave = repmat(ave,numel(amplitude),1);
                        otherwise
                            error('Number of input ave must equal number of input amplitudes, or be equal to 1')
                    end
            end
            varargout = cell(1, nargout);
            
            for i = 1:numel(amplitude)
                A = amplitude(i);   %amplitude of given wave
                switch shape{i} %create fields
                    case 'sin'
                        varargout{i} = A*sin([1:npts]*2*pi/pts_per_cycle) + ave(i);
                    case 'cos'
                        varargout{i} = A*cos([1:npts]*2*pi/pts_per_cycle) + ave(i);
                    case 'ramp'
                        temp = A*(linspace(0,1,pts_per_cycle) - 1/2) + ave(i);
                        for j = 1:ceil(npts/pts_per_cycle)
                            temp = [temp(:) temp(:)];
                        end
                        varargout{i} = temp(1:npts);
                    case 'const'
                        varargout{i} = A*ones(1,npts);
                end
            end
            
        end
        
        %% MH Measurements: Centroid
        function [cg,area]=centroid(x1,y1)
            %CENTROID calculates the geometric centroid of input xy data
            %sets
            %the various input fields.
            %   Inputs:
            
            % make a copy of the pts offset by 1
            x2=x1([2:end 1]);
            y2=y1([2:end 1]);
            
            % compute partial terms
            da = x1.*y2 - x2.*y1;
            dx = (x2 + x1) .* da;
            dy = (y2 + y1) .* da;
            
            % sum
            a = sum(da);
            x = sum(dx);
            y = sum(dy);
            
            % use those to compute signed area & centroid
            area = a / 2;
            cg = [x/(3*a),y/(3*a)];
        end
        
        
        
    end
    
end

