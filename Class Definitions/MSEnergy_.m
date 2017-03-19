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
                obj = Demag_Anisotropy(obj);
                
                %Crystalline anisotropy
                obj = MagnetoCrystalline_Anisotropy(obj);
                
                %magnetoelastic anisotropy (OHandley page 259 (232) )
                obj.U_ME = @(m1,m2,m3,s1,s2,s3,s4,s5,s6) ...
                    obj.MP.B_me(1).*((m1.^2-1/3).*s1 + (m2.^2-1/3).*s2 + (m3.^2-1/3).*s3) + ...
                    obj.MP.B_me(2).*(m2.*m3.*s4 + m1.*m3.*s5 + m1.*m2.*s6);
                
                %uniaxial anisotropy
                obj.U_uni = @(m1,m2,m3) 0*m1;
                
                %exchange bias anisotropy
                obj.U_EB = @(m1,m2,m3) 0*m1;
                
                %PMA energy E=Km3^2
                obj.U_PMA = @(m3) ...
                    obj.MP.Kpma.*(m3.^2);
                
            else
                error('MSEnergy must be constructed with an "MSProperties_" class argument (MSEnergy_(MSProperties_))')
            end
            
            
        end
        %% Total Energy
        function obj = Compute_Total_Energy(obj)                                    
            %Sort out energy terms from property list
            Property_list = properties(obj);
            Energy_ind = ~cellfun(@isempty,strfind(Property_list, 'U_'));
            Energy_list = Property_list(Energy_ind);
            Energy_list = setdiff(Energy_list,'U_total');
            
            %Determine input variables for all energy expressions
            variable_set = cell(1,numel(Energy_list));
            for i = 1:numel(Energy_list)
                temp_func = obj.(Energy_list{i});
                fstr = func2str(temp_func);                
%                 expr = '[\(,]{1}(\w*)[,\)]{1}';
                expr = '[\(,]{1}([a-zA-Z]+[1-6]?)[,\)]{1}';
                [~,endInd] = regexp(fstr,expr );
                [~,tokens1] = regexp(fstr,expr,'match','tokens');                
                [~,tokens2] = regexp(fstr(endInd(1):end),expr,'match','tokens');  
                
                %combine tokens
                temp_vars = {};
                nVars = numel(tokens1) + numel(tokens2);                
                temp_vars(1:2:nVars) = [tokens1{:}] ;
                temp_vars(2:2:nVars) = [tokens2{:}] ;                     
                variable_set{i} = temp_vars;
            end
            
            %Determine unique variables
            unique_variables = unique([variable_set{:}]);
            
            %assemble total energy 
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

%% Function for class only
function obj = Demag_Anisotropy(obj)

N = DemagFactor(obj.MP.Shape,obj.MP.Dims);
obj.DemagFactors = N;
obj.U_Demag = @(m1,m2,m3) 1/2.*obj.MP.mu0*obj.MP.Ms^2*.../
    (N(1)*m1.^2 + N(2)*m2.^2 + N(3)*m3.^2);
end

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
