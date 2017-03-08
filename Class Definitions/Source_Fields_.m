classdef Source_Fields_ <  matlab.mixin.Copyable
    %SOURCE_FIELDS_ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        static
        dynamic
    end
    
    methods
        %% Constructor
        function obj = Source_Fields_(varargin)
            switch isempty(varargin)
                case 1 %default constructor: create empty fields
                    obj.static = Static_Fields_(varargin{1:end-1});
                    obj.dynamic = Dynamic_Fields_(varargin{1:end-1});
                case 0 %assign specified fields
                    switch varargin{end}
                        case 'static' 
                            obj.static = Static_Fields_(varargin{1:end-1});
                        case 'dynamic'
                            obj.dynamic = Dynamic_Fields_(varargin{1:end-1});
                        case 'update'                            
                            obj.static = Static_Fields_(varargin{1}{:} );
                            obj.dynamic = Dynamic_Fields_(varargin{2}{:} );
                    end
            end
        end
        
        
        %% Update
        function obj = Update(obj,particle)
            static_fields = assign_props(particle.SourceFields.static);
            dynamic_fields = assign_props(particle.SourceFields.dynamic);                       
            obj = Source_Fields_(static_fields,dynamic_fields,'update');
        end
        
    end
    
    %% Copy function
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all four properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the objects
            cpObj.static = copy(obj.static);
            cpObj.dynamic = copy(obj.dynamic);
        end
    end
end

%% Local functions
function Name_Value = assign_props(fields)
props = fields.Property_List;
count = 1;
Name_Value = cell(1,numel(props)*2);
for i = 1:numel(props)
    Name_Value{count} = props{i};
    Name_Value{count+1} = fields.(props{i});
    count = count + 2;
end
end
