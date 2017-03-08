classdef MSEffectiveFields_< matlab.mixin.Copyable
    %ENERGY_EXPRESSIONS_' Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(GetAccess = protected)
        %external field anisotropy
        H1_Zeeman
        H2_Zeeman
        H3_Zeeman
        
        %Demag anisotropy
        H1_Demag
        H2_Demag
        H3_Demag
        
        %Crystalline anisotropy
        H1_MCA
        H2_MCA
        H3_MCA
        
        %magnetoelastic anisotropy
        H1_ME
        H2_ME
        H3_ME
        
        %exchange bias anisotropy
        H1_EB
        H2_EB
        H3_EB
        
        %total energy
        H1_total
        H2_total
        H3_total
    end
    
    properties (Hidden = true)
        %Material Properties
        MP
    end
    
    %% Methods
    methods
        %Constructor
        function obj = MSEffectiveFields_(varargin)
            if nargin ==1 && isa(varargin{1},'MSParticle_') %Default Constructor                
                %Assign material properties
                particle = varargin{1};
                obj.MP = particle.Properties;
                
                %% Effective Field Expressions                
                %external field anisotropy
                obj.H1_Zeeman = @(H1,H2,H3) H1;
                obj.H2_Zeeman = @(H1,H2,H3) H2;
                obj.H3_Zeeman = @(H1,H2,H3) H3;
                
                %Demag anisotropy
                N = particle.Energy.DemagFactors;
                obj.H1_Demag = @(m1,m2,m3) -obj.MP.Ms*N(1)*m1;
                obj.H2_Demag = @(m1,m2,m3) -obj.MP.Ms*N(2)*m2;
                obj.H3_Demag = @(m1,m2,m3) -obj.MP.Ms*N(3)*m3;
                
                %Crystalline anisotropy
                obj = MagnetoCrystalline_Anisotropy(obj);
                
                %magnetoelastic anisotropy
                obj.H1_ME = @(m1,m2,m3,s1,s2,s3,s4,s5,s6) -1/(obj.MP.mu0*obj.MP.Ms).*( ...
                    obj.MP.B_me(1).*( (2.*m1).*s1 ) +  obj.MP.B_me(2).*(m3.*s5 + m2.*s6) );
                obj.H2_ME = @(m1,m2,m3,s1,s2,s3,s4,s5,s6) -1/(obj.MP.mu0*obj.MP.Ms).*( ...
                    obj.MP.B_me(1).*( (2.*m2).*s2 ) +  obj.MP.B_me(2).*(m3.*s4 + m1.*s6) );
                obj.H3_ME = @(m1,m2,m3,s1,s2,s3,s4,s5,s6) -1/(obj.MP.mu0*obj.MP.Ms).*( ...
                    obj.MP.B_me(1).*( (2.*m3).*s3 ) +  obj.MP.B_me(2).*(m2.*s4 + m1.*s5) );
                
                %exchange bias anisotropy
                obj.H1_EB = @(m1,m2,m3,H1,H2,H3) 0*m1;
                obj.H2_EB = @(m1,m2,m3,H1,H2,H3) 0*m1;
                obj.H3_EB = @(m1,m2,m3,H1,H2,H3) 0*m1;
                
            else
                error('MSEnergy must be constructed with an "MSProperties_" class argument (MSEnergy_(MSProperties_))')
            end
            
            %compute total effective field
            obj = Compute_Total_Field(obj);
        end
        %% Total Effective Field
        function obj = Compute_Total_Field(obj)
            obj = assign_fields(obj,'H1');
            obj = assign_fields(obj,'H2');
            obj = assign_fields(obj,'H3');
        end
        
        %% Update
        function obj = Update(obj,particle)
            obj = MSEffectiveFields_(particle);
            obj = Compute_Total_Field(obj);
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

%% Function for class only
function obj = MagnetoCrystalline_Anisotropy(obj)
crystal = obj.MP.Crystal;
K = obj.MP.K_mca;
switch crystal
    case 'Cubic'
        obj.H1_MCA = @(m1,m2,m3) -1/(obj.MP.mu0*obj.MP.Ms).* (K(1).*(2*m1 .* m2.^2 + 2.* m1 .* m3.^2 ) +...
            K(2).*(2 .* m1 .* m2.^2 .* m3.^2) );
        obj.H2_MCA = @(m1,m2,m3) -1/(obj.MP.mu0*obj.MP.Ms).* (K(1).*(2*m2 .* m1.^2 + 2.* m2 .* m3.^2 ) +...
            K(2).*(2 .* m2 .* m1.^2 .* m3.^2) );
        obj.H3_MCA = @(m1,m2,m3) -1/(obj.MP.mu0*obj.MP.Ms).* (K(1).*(2*m3 .* m1.^2 + 2.* m3 .* m2.^2 ) +...
            K(2).*(2 .* m3 .* m1.^2 .* m2.^2) );
    case 'Hexagonal'
        sTh2 = @(m3) 1-m3.^2;    %sin(theta)^2
        phi = @(m1,m2) atan2(m2,m1);
        obj.H1_MCA = @(m1,m2,m3) -1/(obj.MP.mu0*obj.MP.Ms).* ( K(3).*sTh2(m3).^3.* 6 .* sin(6.*phi(m1,m2))*(-m1./sqrt(m1.^2 + m2.^2 +eps)) );
        obj.H1_MCA = @(m1,m2,m3) -1/(obj.MP.mu0*obj.MP.Ms).* ( K(3).*sTh2(m3).^3.* 6 .* sin(6.*phi(m1,m2))*( m2./sqrt(m1.^2 + m2.^2 +eps)) );
        obj.H3_MCA = @(m1,m2,m3) -1/(obj.MP.mu0*obj.MP.Ms).* ( -2*K(1).* m3 + 2*K(2).*sTh2(m3).*(-2.*m3) + ...
            3.*K(3).*sTh2(m3).^2.*(-2.*m3).*cos(6.*phi(m1,m2)) );
    case 'Amorphous'
        obj.H1_MCA = @(m1,m2,m3) 0*(m1+m2+m3);   %leave values so it returns the right size
        obj.H2_MCA = @(m1,m2,m3) 0*(m1+m2+m3);
        obj.H3_MCA = @(m1,m2,m3) 0*(m1+m2+m3);
    case 'Poly'
        obj.H1_MCA = @(m1,m2,m3) 0*(m1+m2+m3);   %leave values so it returns the right size
        obj.H2_MCA = @(m1,m2,m3) 0*(m1+m2+m3);
        obj.H3_MCA = @(m1,m2,m3) 0*(m1+m2+m3);
    otherwise
        error('%s is not an accepted crystal type. Allowable crystals are: \n%s',crystal,strjoin(obj.MP.Allowable_Crystals,', ' ))
        
end


end

function obj = assign_fields(obj,name)
%Sort out energy terms from property list
Property_list = properties(obj);
ind = ~cellfun(@isempty,strfind(Property_list, [name,'_']));
list = Property_list(ind);
list = setdiff(list,[name,'_total']);

%Determine input variables for all energy expressions
variable_set = cell(1,numel(list));
for i = 1:numel(list)
    temp_func = obj.(list{i});
    fstr = func2str(temp_func);
    %expr = '[\(,]{1}(\w*)[,\)]{1}';
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
eval_str = ['obj.',list{1},['(',strjoin(variable_set{1},','),')']];
for i = 2:numel(list)
    eval_str = [eval_str, ' + ',['obj.',list{i},['(',strjoin(variable_set{i},','),')']] ];
end
eval_str = ['obj.',name,'_total = @(',strjoin(unique_variables,','),') ', eval_str,';'];
eval(eval_str);
end