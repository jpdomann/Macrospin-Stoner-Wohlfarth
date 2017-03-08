classdef MSState_ < matlab.mixin.Copyable
    %MSSTATE_ Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Properties
    properties (SetAccess = protected, GetAccess = public)
        m1          %initial orientation
        m2
        m3
        
        Stat_m1     %Static Data
        Stat_m2
        Stat_m3
        Stat_theta
        Stat_phi
        
        Dyn_m1      %Dynamic Data
        Dyn_m2
        Dyn_m3
        Dyn_theta
        Dyn_phi
    end
    properties (Access = public, Constant = true)
        Property_List = {'m','m1','m2','m3'};
    end
    
    %% Methods
    methods
        %default constructor
        function obj = MSState_(varargin)
            switch nargin <= 1
                case 1 %default constructor
                    obj.m1 = 1;
                    obj.m2 = 0;
                    obj.m3 = 0;
                case 0
                    %particle is first passed object
                    particle = varargin{1};
                    
                    %separate input into property / value pairs
                    %these should only contain initial magnetization
                    %directions
                    [prop, val] = PropertyValue(varargin(2:end));
                    %add Init_ if prop's are of form 'm#'
                    
                    %Check for wrong inputs
                    Allowable_Props = obj.Property_List;
                    NotAllowable = ~ismember(prop,Allowable_Props);
                    if any(NotAllowable)
                        error_msg1 = sprintf('"%s" is not a valid property \nValid properties are:\n', prop{NotAllowable} );
                        error_msg2 = sprintf('%s\n',strjoin(Allowable_Props,', ' ));
                        error('%s%s',error_msg1,error_msg2)
                    end
                    
                    %Assign prescribed properties
                    for i = 1:numel(val)
                        switch strcmp(prop{i},'m')
                            case 0
                                obj.(prop{i}) = val{i};
                            case 1
                                if numel(val{i})~= 3 
                                    error('When changing ''m'' all at once, the input data must be a 1x3 or 3x1 vector'); 
                                end
                                obj.m1 = val{i}(1);
                                obj.m2 = val{i}(2);
                                obj.m3 = val{i}(3);
                        end
                    end
                    
                    %Ensure initial state is consistent (abs(m)=1)
                    obj = state_check(obj);
                    
                    %                     %initialize static orientation
                    %                     obj.Stat_m1 = obj.Init_m1 * ones(size(particle.SourceFields.static.H1));
                    %                     obj.Stat_m2 = obj.Init_m2 * ones(size(particle.SourceFields.static.H1));
                    %                     obj.Stat_m3 = obj.Init_m3 * ones(size(particle.SourceFields.static.H1));
                    %
                    %                     %Assign dynamic properties
                    %                     obj.Dyn_m1 = obj.Init_m1 * ones(size(particle.SourceFields.dynamic.t));
                    %                     obj.Dyn_m2 = obj.Init_m2 * ones(size(particle.SourceFields.dynamic.t));
                    %                     obj.Dyn_m3 = obj.Init_m3 * ones(size(particle.SourceFields.dynamic.t));
            end
            
        end
        
        %% Set function
        %assign values to the stationary and dynamic m values
        function obj = set_state(obj,dat,name)
            obj.(name) = dat;
        end
        
        
        %% Update function
        %         function obj = Update(obj,particle)
        %             obj = MSState_(particle.Properties );
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

%% Local functions
function obj = state_check(obj)
%Input direction cosines
for i = 1:3
    switch isempty( eval(['obj.m',num2str(i)]) )
        case 1
            a(i) = 0;
        case 0
            a(i) = eval(['obj.m',num2str(i)]);
    end
end

switch all(a==0)
    case 0  %assign magnetization
        m = a ./ norm(a);
    case 1 %assign a default orientation
        m = [1 0 0];
end
%output state
obj.m1 = m(1);
obj.m2 = m(2);
obj.m3 = m(3);

end