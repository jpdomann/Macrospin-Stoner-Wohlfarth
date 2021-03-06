classdef MS < handle
    %MS Is a utility class that holds useful functions to be accessed
    %by multiple other classes
    
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
        
        %% MH Measurements: Centroid / Heb
        function [cg,area]=centroid(x1,y1)
            %CENTROID calculates the geometric centroid of input xy data
            %sets
            %the various input fields.
            %   Inputs: x1, and y1 are the input x/y data
            
            %input check
            if all(size(x1) ~= size(y1)); error('input x and y must be the same size'); end
            
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
        
        function [ Hc,Heb, Mc, Meb ] = HebANDHc( H,M )
            %HEBANDHC Summary of this function goes here
            %   Detailed explanation goes here
            
            %Break into rising and falling segments
            I_min = find(H == min(H));
            I_min = I_min(1);
            
            if I_min == 1
                I_min = find(H == max(H));
                I_min = I_min(1);
                Hrise = H(1:I_min);
                Mrise = M(1:I_min);
                Hfall = H(I_min:end);
                Mfall = M(I_min:end);
            else
                Hfall = H(1:I_min);
                Mfall = M(1:I_min);
                Hrise = H(I_min:end);
                Mrise = M(I_min:end);
            end
            
            %Ensure monotonic behavior of data (add / subtract eps
            Hfall = MS.makeMonotonic(Hfall);
            Mfall = MS.makeMonotonic(Mfall);
            Hrise = MS.makeMonotonic(Hrise);
            Mrise = MS.makeMonotonic(Mrise);
            
            %Fit an interpolation curve for each region
            fit = 'linear';
            n=1e4;
            hfall = linspace(Hfall(1),Hfall(end),n);
            mfall = interp1(Hfall,Mfall,hfall,fit,'extrap');
            
            hrise = linspace(Hrise(1),Hrise(end),n);
            mrise = interp1(Hrise,Mrise,hrise,fit,'extrap');
            
            %Find crossing locations
            Ifall = find(mfall<0,1);
            Irise = find(mrise>0,1);
            
            %Coercive field
            Hc(1) = hfall(Ifall);
            Hc(2) = hrise(Irise);
            Mc(1) = mfall(Ifall);
            Mc(2) = mrise(Irise);
            %Exchange Bias
            Heb = mean(Hc);
            Meb = mean(Mc);
            
            % %plot
            % plot(Hfall,Mfall,'r',Hrise,Mrise,'b')
            % hold on
            % plot(hfall,mfall,'b.',hrise,mrise,'r.')
            % plot(Hc,Mc,'g*',Heb,Meb,'g*')
            % grid on
            
        end        

        %% Check Numeric Arrays
        function [ obj ] = Check_Numeric_Arrays(obj,Property_List,err_name)
            %CHECK_NUMERIC_ARRAYS provides a check on the values of numeric
            %data arrays.
            %   Loop over all the arrays in a given property list
            
            %ensure all values are numeric and not NaN
            sz = zeros(1,numel(Property_List));
            for i = 1:numel(Property_List)
                val = obj.(Property_List{i});
                sz(i) = numel(val);
                nanFlag = isnan(val);
                numFlag = ~isnumeric(val);
                if any(nanFlag)
                    error('NaN numbers encountered in %s',Property_List{i})
                elseif any(numFlag)
                    error('Non-numbers encountered in %s',Property_List{i})
                end
            end
            %make sure all arrays are the same size
            switch all(sz == sz(1));
                case 0
                    error('Not all %s are the same size',err_name)
                case 1
                    
            end
            
        end
        
        %% Demag Factors
        function N = DemagFactor(shape,Dims)
            %DEMAGFACTOR calculates the demagnetization factors for ellipsoidal and
            %rectangular shapes.
            %   To do: enable the user to rotate the ellipse with a set of Euler angles
            
            [Dims_sorted, Sort_ind] = sort(Dims,'descend');
            
            a = Dims_sorted(1);
            b = Dims_sorted(2);
            c = Dims_sorted(3);
            if (a <0 || b<0 || c<0)
                warning('A negative dimension was entered')
                N = [0 0 0];
                return
            end
            
            switch shape
                case 'Ellipse'
                    if (a == b) && (b == c) %sphere
                        N1 = 1/3;
                        N2 = 1/3;
                        N3 = 1/3;
                    elseif (a==b) && a>c     %oblate spheroid
                        e = sqrt((a^2-c^2)/(a^2));
                        N3 =1/e^2 .* (1 - sqrt(1-e^2)/e *asin(e));
                        N1 = 0.5 .* (1-N3);
                        N2 = N1;
                    elseif (a == b) && (a < c)  %prolate spheroid
                        e = sqrt((c^2-a^2)/(c^2));
                        N3 =(1-e^2)/e^2 .* (1/(2*e).*log((1+e)/(1-e)) -1 );
                        N1 = 0.5 .* (1-N3);
                        N2 = N1;
                    elseif (a > b) && (b == c)  %prolate spheroid
                        m = a/c;
                        n = m^2-1;
                        N1 = 1/n * ( m /(2*sqrt(n)) * log((m+sqrt(n))/(m-sqrt(n))) - 1 );
                        N2 = m/(2*n) * (m - 1/(2*sqrt(n))* log((m+sqrt(n))/(m-sqrt(n))) );
                        N3 = N2;
                    elseif (a >= b) && (b >= c)
                        %"Depolartization Tensor Method - 2nd Ed."
                        B = b/a;    B2 = 1-B^2;
                        G = c/a;    G2 = 1-G^2;
                        k = sqrt(B2 / G2 )^2;
                        phi = acos(G);
                        N1 = B*G /( B2*sqrt(G2) )*( ellipticF(phi,k) - ellipticE(phi,k) );
                        N2 = - G^2/(B^2 - G^2) +...
                            B*G*sqrt(G2) / ( B2*(B^2-G^2) ) * ellipticE(phi,k) -...
                            B*G/( B2*sqrt(G2) ) * ellipticF(phi,k);
                        N3 = B^2/(B^2-G^2) -...
                            B*G/( (B^2-G^2)*sqrt(G2) ) * ellipticE(phi,k);
                        
                        %             phi = asin((a^2-c^2)/a^2);
                        %             k = sqrt((a^2 - b^2) / (a^2 - c^2))^2;
                        %             N1 = a*b*c / (sqrt(a^2 - c^2) * (a^2 - b^2) ) * (-ellipticE(phi,k) + ellipticF(phi,k));
                        %             N2 = -a*b*c^2 / ( a*b*(b^2-c^2) ) ...
                        %                  +a*b*c*sqrt(a^2-c^2) / ( (a^2-b^2)*(b^2-c^2) )* ellipticE(phi,k) ...
                        %                  -a*b*c / ( sqrt(a^2-c^2)*(a^2-b^2) ) * ellipticF(phi,k);
                        %             N3 =  a*b^2*c/(a*c*(b^2-c^2)) ...
                        %                  -a*b*c/(sqrt(a^2-c^2)*(b^2-c^2)) * ellipticE(phi,k) ;
                    else
                        error('Demag error: For a general ellipsoid, dimensions must be ordered [a,b,c], where a>=b>=c')
                        
                    end
                    
                case 'Rectangle'
                    m = @(a,b,c) sqrt(a^2 + b^2 + c^2);
                    D = @(a,b,c) ...
                        (b^2-c^2)/(2*b*c)*log( (m(a,b,c) - a) / ( m(a,b,c) + a ) ) + ...
                        (a^2-c^2)/(2*a*c)*log( (m(a,b,c) - b) / ( m(a,b,c) + b ) ) + ...
                        b/(2*c)*log( (m(a,b,0) +a )/( m(a,b,0)-a ) ) + ...
                        a/(2*c)*log( (m(a,b,0) +b )/( m(a,b,0)-b ) ) + ...
                        c/(2*a)*log( (m(0,b,c) -b )/( m(0,b,c)+b ) ) + ...
                        c/(2*b)*log( (m(a,0,c) -a )/( m(a,0,c)+a ) ) + ...
                        2*atan( (a*b) / ( c * m(a,b,c) ) )+ ...
                        (a^3 + b^3 -2*c^3)/(3*a*b*c) + ...
                        (a^2 + b^2 -2*c^2)/(3*a*b*c)*m(a,b,c) + ...
                        c/(a*b)*( m(a,0,c) + m(0,b,c) ) - ...
                        (m(a,b,0)^3  +m(0,b,c)^3 + m(a,0,c)^3)/(3*a*b*c);
                    
                    N3 = D(a,b,c)/(pi);
                    N1 = D(b,c,a)/(pi);
                    N2 = D(c,a,b)/(pi);
                otherwise
                    shapes = 'Ellipse, Rectangle';
                    error(sprintf('%s is not an allowable shape. \nAllowable shapes are: %s',shape,shapes))
            end
            
            check = N1+N2+N3 -1 < 1e-6;
            if ~check
                error('Demag factors don''t add to 1!')
            end
            
            %return to original order
            N(Sort_ind(1))= N1;
            N(Sort_ind(2))= N2;
            N(Sort_ind(3))= N3;
            
        end
        
        %% Demag Anisotropy
        function obj = Demag_Anisotropy(obj)
            
            N = MS.DemagFactor(obj.MP.Shape,obj.MP.Dims);
            obj.DemagFactors = N;
            obj.U_Demag = @(m1,m2,m3) 1/2.*obj.MP.mu0*obj.MP.Ms^2*.../
                (N(1)*m1.^2 + N(2)*m2.^2 + N(3)*m3.^2);
        end
        
        %% Magnetocrystalline Anisotropy        
        function obj = MagnetoCrystalline_Anisotropy(obj)
            crystal = obj.MP.Crystal;
            K = obj.MP.K_mca;
            switch crystal
                case 'Cubic'
                    obj.U_MCA = @(m1,m2,m3) K(1).*(m1.^2 .* m2.^2 + m1.^2 .* m3.^2 + m2.^2 .* m3.^2 ) +...
                        K(2).*(m1.^2 .* m2.^2 .* m3.^2);     %Cullity - pg 215
                case 'Hexagonal'
                    sTh2 = @(m3) 1-m3.^2;    %sin(theta)^2
                    phi = @(m1,m2) atan2(m2,m1);
                    obj.U_MCA = @(m1,m2,m3) K(1).*sTh2(m3) + K(2).*sTh2(m3).^2 + ...
                        K(3).*sTh2(m3).^3.*cos(6.*phi(m1,m2));
                case 'Amorphous'
                    obj.U_MCA = @(m1,m2,m3) 0*(m1+m2+m3);   %leave values so it returns the right size
                case 'Poly'
                    obj.U_MCA = @(m1,m2,m3) 0*(m1+m2+m3);   %leave values so it returns the right size
                otherwise
                    error('%s is not an accepted crystal type. Allowable crystals are: \n%s',crystal,strjoin(obj.MP.Allowable_Crystals,', ' ))                    
            end                        
        end

        %% Anonymous Function Argument / Variable List
        function [args,list] = Energy_argument_list(obj)
            %ENERGY_ARGUMENT_LIST takes an input MSEnergy_ object, and returns the
            %argument list for each of the energy terms
            
            %Sort out energy terms from property list
            Property_list = properties(obj);
            if any(cell2mat(strfind(Property_list,'U_')))
                ind = ~cellfun(@isempty,strfind(Property_list, 'U_'));
                list = Property_list(ind);
                list = setdiff(list,'U_total');
            elseif any(cell2mat(strfind(Property_list,'H1')))
                ind = ~cellfun(@isempty,strfind(Property_list, '_total'));
                list = Property_list(~ind);
            end
            %Determine input variables for all energy expressions
            args = cell(1,numel(list));
            for i = 1:numel(list)
                args{i} = MS.argument_list(obj.(list{i}));
            end
        end
        
        function args = argument_list(func)
            %ARGUMENT_LIST takes an input function handle, and returns the
            %arugment list        
            switch isnumeric(func)
                case 1
                    args = [];
                case 0
                    args = strsplit(regexp(func2str(func), '(?<=^@\()[^\)]*', 'match', 'once'), ',');
            end
        end
        
        %% Utility Function - Make Monotonic
        function [ x_mono ] = makeMonotonic( x )
            %MAKEMONOTONIC makes sure data is either monotonically increasing or
            %decreasing
            
            
            min_val = 1e-13;
            difference = diff(x);
            switch x(1) < x(end)
                case 1 % data is increasing
                    if any(difference <=0)
                        I = difference <= 0;
                        I = [I(:);false];
                        if find(I == 1) > 1
                            x(I) = x(I)+min_val;
                        elseif find(I == 1) == 1
                            x(I) = x(I)-min_val;
                        end
                    end
                case 0 %data is decreasing
                    if any(difference >=0)
                        I = difference >= 0;
                        I = [I(:);false];
                        if find(I == 1) > 1
                            x(I) = x(I)-min_val;
                        elseif find(I == 1) == 1
                            x(I) = x(I)+min_val;
                        end
                    end
            end
            
            x_mono = x;
            
        end
        
        %%
        
        %%
               
    end
    
end

