classdef MSProperties_ < matlab.mixin.Copyable
    %MSPROPERTIES Contains the properties for a MSParticle
    %   Details: It is inherited from the handle class, enabling calls to
    %   the MSParticle to change it's properties
    
    properties
        Location        %XYZ location of particle center
        Shape           %Type of shape: ellipse or rectangle
        Dims            %Dimensions rectangle:(l,w,h), ellipse:(a,b,c)
        Volume
        Mat_Name        %Material Name
        Ms              %Saturation magnetization (A/m = emu/cc *1000)
        Crystal         %Type of crystal: cubic, hex, amorphous
        K_mca           %MCA coefficients
        B_me            %Magnetoelastic Coefficients
        Alpha           %Gilbert Damping 
        Kuni            %uniaxial anisotropy coefficient
        dir_uni         %uniaxial anisotropy axis
        Keb             %exchange bias anisotropy coefficient
        dir_eb          %exchange bias direction
        Ks              %surface PMA coefficient
        Kpma            %PMA coefficient Kpma=Ks/t, t is thickness
    end
    properties (Constant)
        mu0 = 4*pi*1e-7;    %[H/m] vacuum permeability
        Gamma_e =1.760859e11  %[1/(s*T)] Electron Gyromagnetic ratio        
        
        Property_List = {...%Property Descriptions
            'Property', 'Description',                          'Units';...
            'Location'  'XYZ location of particle''s center'    'm';...
            'Shape',    'Material Shape',                       '';...
            'Dims',     'Dimensions [a,b,c]',                   'm';...  %2a,2b,2c for an ellipse (outer dimensions, not semi-axes)
            'Mat_Name', 'Material Name',                        '';...
            'Ms',       'Saturation Magnetization',             'A/m';...
            'Crystal',  'Crystal type (cubic,hex,amporphous)',  '';
            'K_mca',    'MCA Coefficients [K1,K2,...]',         'J/m^3';...
            'B_me',     'Magnetoelastic Coefficients',          'J/m^3';...
            'Alpha',    'Gilbert Damping',                      '';...
            'Kuni',     'Uniaxial anisotropy coefficient',      'J/m^3';...
            'dir_uni',  'Uniaxial anisotropy axis',             '';...
            'Keb',      'Exchange bias coefficient',            'J/m^3';...
            'dir_eb',   'Exchange bias direction',              '';...
            'Ks',       'surface coefficient',                  '';...
            'Kpma',     'PMA coefficient',                      '';...
            };
        Available_Mats = 'Iron, Nickel, Cobalt, CoFeB';
        Allowable_Shapes = {'Ellipse','Rectangle'};
        Allowable_Crystals = {'Cubic','Hexagonal','Amorphous','Poly'};
    end
    properties (SetAccess = protected)
        Gamma       %= mu0*Gamma_e;
        Factor_1    % = Gamma / (1 + Alpha^2);
        Factor_2    % = Gamma * Alpha / (1 + Alpha^2);
    end
    methods
        %Default Constructor
        function obj = MSProperties_(varargin)
            if nargin <= 1      %load properties for specified material                                
                %Parse input
                if nargin == 0
                    input = 'Nickel'; %Nickel is default material
                else
                    input = varargin{1};
                end
                
                %Check for string input
                if ischar(input) && any(strfind(obj.Available_Mats,input))
                    obj.Mat_Name = input;
                else
                    error('Single input needs to be valid material name:\n%s',obj.Available_Mats')
                end
                
                %Assign pre-defined properties
                obj = DefaultProperties(obj,input);
                
            else
                %separate input into property / value pairs
                [prop, val] = MS.PropertyValue(varargin);
                
                %Check for wrong inputs
                Allowable_Props = {obj.Property_List{2:end,1}};
                NotAllowable = ~ismember(prop,Allowable_Props);
                if any(NotAllowable)
                    error_msg1 = sprintf('"%s" is not a valid property \nValid properties are:\n', prop{NotAllowable} );
                    error_msg2 = sprintf('%s\n',strjoin(Allowable_Props,', ' ));
                    error('%s%s',error_msg1,error_msg2)
                end
                
                %If a material name was provided, load default values if
                %the material already exists
                Ind = strcmp(prop,'Mat_Name');
                obj = DefaultProperties(obj,val{Ind});
                
                
                %Assign prescribed properties
                for i = 1:numel(val)
                    obj.(prop{i}) = val{i};
                end
                
                %Check MSParticle is in a valid state
                obj = check_Properties(obj);
                                
            end
            obj.Gamma = obj.mu0*obj.Gamma_e;
            obj.Factor_1 = obj.Gamma / (1 + obj.Alpha^2);
            obj.Factor_2 = obj.Gamma * obj.Alpha / (1 + obj.Alpha^2);
        end
        
        function obj = check_Properties(obj)
            %Check MSParticle is in a valid state
            obj = check_Geom(obj);
            obj = check_Ms(obj);
            obj = check_Crystal(obj);
            obj = check_Factors(obj);
            obj = check_PMA(obj);
            obj = check_K(obj);
            obj = check_Directions(obj);
        end

        function obj = Load_Default_Material(obj,MatName)
            %verify proper input
            if nargin < 2
                error('Not enough inputs. Enter material name:\n%s',obj.Available_Mats')
            elseif nargin > 2
                error('Too many inputs. Enter single material name:\n%s',obj.Available_Mats')
            end
            
            %load a set of default material properties
            obj = DefaultProperties(obj,MatName);
            
        end
    end
    
    %% Copy function
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all four properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);

        end
    end
    
end

%% Local methods: Available to Class only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = DefaultProperties(obj,input)

%set remaining properties (these will be overwritten later if the user
%specified them)
obj.Location = [0,0,0];
obj.Shape = 'Ellipse';
obj.Dims = [1 1 1];

switch input
    case 'Nickel'
        obj.Mat_Name = 'Nickel';
        obj.Ms = 484*1e3;
        obj.Crystal = 'Cubic';
        obj.K_mca = [-4.5E3,-2.3E3];	%OHandley page 216 (189)
        obj.B_me = [6.2e6 4.3e6];       %OHandley page 257 (230)
        obj.Alpha = 0.01;
        obj.Ks = -2.6e-4;        
    case 'Iron'
        obj.Mat_Name = 'Iron';
        obj.Ms = 1714*1e3;
        obj.Crystal = 'Cubic';
        obj.K_mca = [48E3,-10E3];	%OHandley page 216 (189)
        obj.B_me = [-2.9e6 2.9e6];  %OHandley page 257 (230)
        obj.Alpha = 0.01;
        obj.Ks = 0;
    case 'Cobalt'
        obj.Mat_Name = 'Cobalt';
        obj.Ms = 1422*1e3;
        obj.Crystal = 'Hexagonal';
        obj.K_mca = [450E3,150E3, 0];	%Cullity page 241(227)
        obj.B_me = [6e6 13e6];          %OHandley page 257 (230)
        obj.Alpha = 0.01;
        obj.Ks = 0;
    case 'CoFeB'
        obj.Mat_Name = 'CoFeB';
        obj.Ms = 1.2e6; %Chen - 1.11e6, QC 1.2e6;
        obj.Crystal = 'Poly';
        obj.K_mca = [0,0];              %negligible for poly-crystal
        obj.B_me = [-2.7e7 -2.7e7];     %use lambda_s = 150ppm
        obj.Alpha = 0.01;
        obj.Ks = -1.5e-3;
    case 'Galfenol'
        obj.Mat_Name = 'Galfenol';
        obj.Ms = 1.456e6;                 %
        obj.Crystal = 'Poly';
        obj.K_mca = [0,0];              %negligible for poly-crystal
        obj.B_me = [-2.1635e+7 -4.3269e+7];     %use lambda_s = 150ppm
        obj.Alpha = 0.01;
        obj.Ks = 0e-3;
    case 'Galfenol-single-crystal'
        obj.Mat_Name = 'Galfenol-single-crystal';
        obj.Ms = 1.456e6;                 %
        obj.Crystal = 'Cubic';
        obj.K_mca = [13e3, -90e3];              %negligible for poly-crystal
        obj.B_me = [-9.84e+6 4.797e+6];     %use lambda_s = 150ppm
        obj.Alpha = 0.01;
        obj.Ks = 0e-3;
    otherwise
        error('Single input needs to be valid material name:\n%s',obj.Available_Mats')
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = check_Geom(obj)
%check location
loc = obj.Location;
if iscell(loc) || ~isnumeric(loc)
   error('Location must be a numeric vector, cell arrays are not allowed ') 
end

%check shape
allowable = obj.Allowable_Shapes;
current_val = obj.Shape;
if isempty(current_val)
    error('Must enter a Shape. Allowable shapes are: \n%s',strjoin(allowable,', '))
end
check = ismember(current_val,allowable);
if ~check
    error('%s is not an allowable Shape. Allowable shapes are: \n%s',current_val,strjoin(allowable,', '))
end

%check dimensions
current_dims = obj.Dims;
if numel(current_dims)<3 || any(current_dims <= eps)
    error('When Shape is adjusted, you must enter 3 positive dimensions (Dims) [x,y,z] > %e',eps)
end

%update volume
switch obj.Shape
    case 'Ellipse'
        obj.Volume = 4/3 * pi * obj.Dims(1) * obj.Dims(2) * obj.Dims(3) / 8;
    case 'Rectangle'
        obj.Volume = obj.Dims(1) * obj.Dims(2) * obj.Dims(3);
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = check_Ms(obj)
Ms = obj.Ms;
if isempty(Ms)
    error('Must enter an Ms value. Ms must be positive')
end

check = Ms <= 0;
if check
    error('Ms must be a positive value')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = check_Factors(obj)
Alpha = obj.Alpha;
check = Alpha <= 0;
if check
    error('Alpha must be a positive value')
end
obj.Factor_1 = obj.Gamma / (1 + obj.Alpha^2);
obj.Factor_2 = obj.Gamma * obj.Alpha / (1 + obj.Alpha^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
function obj = check_Crystal(obj)
%check crystal type
allowable = obj.Allowable_Crystals;
current_val = obj.Crystal;
if isempty(current_val)
    error('Must enter a Crystal type. Allowable crystals are: \n%s',strjoin(allowable,', '))
end

check = ismember(current_val,allowable);
if ~check
    error('%s is not an allowable Crystal. Allowable crystals are: \n%s',current_val,strjoin(allowable,', '))
end

%check K values
current_K = obj.K_mca;
%cubic and hex crystals must have K_mca values entered
if isempty(current_K) && any(strcmp(obj.Crystal,{'Cubic','Hexagonal'}))
    error('You must enter K_mca values along with the Crystal type')
end
numK = numel(current_K);
switch obj.Crystal
    case 'Cubic'
        check = numK~=2;
        if check
            error('Cubic anisotropy requires 2 K_mca values [K1 K2]')
        end
    case 'Hexagonal'
        check = numK~=3;
        if check
            error('Hexagonal anisotropy requires 3 K_mca values [K1 K2 K3]')
        end
    case 'Amorphous'
        %auto set K_mca to zero
        obj.K_mca = 0;
        %isotropic magnetoelasticity
        obj.B_me(2) = 2*obj.B_me(1);
    case 'Poly'
        %auto set K_mca to zero
        obj.K_mca = 0;
        %isotropic magnetoelasticity
        obj.B_me(2) = 2*obj.B_me(1);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
function obj = check_PMA(obj)
%update PMA terms
obj.Kpma = obj.Ks/obj.Dims(3);

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = check_K(obj)
%uniaxial anisotropy
switch isempty(obj.Kuni)        
    case 1
        obj.Kuni = 0;
end
%exchange bias anisotropy
switch isempty(obj.Keb)        
    case 1
        obj.Keb = 0;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function obj = check_Directions(obj)
%uniaxial anisotropy directions
obj.dir_uni = anisotropy_direction_check(obj.dir_uni);
%exchange bias anisotropy directions
obj.dir_eb = anisotropy_direction_check(obj.dir_eb);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function anis_direc = anisotropy_direction_check(anis_direc)

switch isempty(anis_direc)
    case 1 %if no direction is provided, set to zero
        anis_direc = zeros(3,1);
    case 0 %if direction is provided, normalize to one
        %check for correct number of dimensions
        if numel(anis_direc) ~= 3; error('Must input a 3 dimensional vector for dir_uni');end 
        %normalize to one
        switch all(anis_direc==0) 
            case 0 %if a non-zero vector is entered
                anis_direc  = anis_direc  ./ norm(anis_direc );
            case 1 %if a zero vector is entered
              anis_direc = zeros(3,1);  
        end
end
%make column vector
anis_direc  = anis_direc(:);

end

