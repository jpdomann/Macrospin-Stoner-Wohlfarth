classdef MSEnergy_< matlab.mixin.Copyable
    %ENERGY_EXPRESSIONS_' Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(GetAccess = protected)        
        U_Zeeman	%external field anisotropy
        U_Demag     %Demag anisotropy
        DemagFactors                
        U_MCA       %Crystalline anisotropy        
        U_ME        %magnetoelastic anisotropy        
        U_EB        %exchange bias anisotropy        
        U_uni       %uniaxial anisotropy        
        U_PMA       %PMA energy (Uniaxial in z-direction)        
        U_total     %total energy
    end
    
    properties (Hidden = true)
        %Material Properties
        MP
    end
    
    %% Methods
    methods
        %Constructor
        function obj = MSEnergy_(varargin)
            if nargin ==1 && strcmp(class(varargin{1}),'MSProperties_')%Default Constructor
                
                %Assign material properties
                obj.MP = varargin{1};
                
                %% Energy Expressions
                
                %external field anisotropy
                obj.U_Zeeman = @(m1,m2,m3,H1,H2,H3) ...
                    -obj.MP.mu0.*obj.MP.Ms.*(m1.*H1+ m2.*H2 + m3.*H3);
                
                %Demag anisotropy
                obj = MS.Demag_Anisotropy(obj);
                
                %Crystalline anisotropy
                obj = MS.MagnetoCrystalline_Anisotropy(obj);
                
                %magnetoelastic anisotropy (OHandley page 259 (232) )
                obj.U_ME = @(m1,m2,m3,s1,s2,s3,s4,s5,s6) ...
                    obj.MP.B_me(1).*((m1.^2-1/3).*s1 + (m2.^2-1/3).*s2 + (m3.^2-1/3).*s3) + ...
                    obj.MP.B_me(2).*(m2.*m3.*s4 + m1.*m3.*s5 + m1.*m2.*s6);
                
                %uniaxial anisotropy dir_uni
                obj.U_uni = @(m1,m2,m3) -obj.MP.Kuni.*(obj.MP.dir_uni(1).*m1+...
                    obj.MP.dir_uni(2).*m2+...
                    obj.MP.dir_uni(3).*m3 ).^2;
                
                %exchange bias anisotropy
                switch isnumeric(obj.MP.Keb) && isempty(MS.argument_list(obj.MP.Keb))
                    case 1
                        obj.U_EB = @(m1,m2,m3) obj.MP.Keb.*(obj.MP.dir_eb(1).*m1+...
                            obj.MP.dir_eb(2).*m2+...
                            obj.MP.dir_eb(3).*m3);
                    case 0 %allow Keb to be a function of other parameters
                        arg_list = MS.argument_list(obj.MP.Keb);
                        str = ['obj.U_EB = @(m1,m2,m3,',strjoin(arg_list,','),') obj.MP.Keb(',strjoin(arg_list,','),').*(obj.MP.dir_eb(1).*m1+ obj.MP.dir_eb(2).*m2+ obj.MP.dir_eb(3).*m3);'];
                        eval(str)
                        
                end
                %PMA energy E=Km3^2
                obj.U_PMA = @(m3) ...
                    obj.MP.Kpma.*(m3.^2);
                
            else
                error('MSEnergy must be constructed with an "MSProperties_" class argument (MSEnergy_(MSProperties_))')
            end
            
            
        end
        %% Total Energy
        function obj = Compute_Total_Energy(obj)       
                                         
            %Determine unique variables
            [variable_set, Energy_list] = MS.Energy_argument_list(obj);     
            unique_variables = unique([variable_set{:}]);
            
            %assemble total energy expression
            eval_str = ['obj.',Energy_list{1},['(',strjoin(variable_set{1},','),')']];
            for i = 2:numel(Energy_list)                
                eval_str = [eval_str, ' + ',['obj.',Energy_list{i},['(',strjoin(variable_set{i},','),')']] ];
            end
            eval_str = ['obj.U_total = @(',strjoin(unique_variables,','),') ', eval_str,';'];
            eval(eval_str);
        end
        
        %% Update
        function obj = Update(obj,particle)
            obj = MSEnergy_(particle.Properties);
            obj = Compute_Total_Energy(obj);
        end
        
%         function obj = Add_Energy_Expression(custom_energy)
%             %have user create an instance of the Custom_MSEnergy_ class.
%             Then feed that class into this function to add additional
%             energy terms
%               custom_energy (name, constants, function(mag, strain) )
%         end                        
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


