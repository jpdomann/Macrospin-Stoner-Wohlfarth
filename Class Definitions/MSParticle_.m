classdef MSParticle_<matlab.mixin.Copyable %copyable handle class
    %MSPARTICLE provides a Stoner-Wohlfarth (macrospin) magnetic particle.
    %The available methods enable the output of quasi-static MH loops for
    %due to multiple anisotropies. Additionally, dynamic magnetization
    %changes can be computed using the Landau-Lifshitz-Gilbert equation of
    %motion (LLG)
    %   Properties (State):
    %   Methods:
    
    %% Propterties
    
    properties %access is public
        Properties      %store relevant properties (material and geometric)
    end
    properties (GetAccess = public, SetAccess = private)
        Energy          %store energy expressions        
        SourceFields    %user supplied stimuli
        State           %state of particle (static / dynamic)
        EffectiveFields %effective magnetic fields for energy expressions
    end
    
    %% Methods
    methods        
        %% Constructors
        function obj = MSParticle_(varargin)
            %Load default values
            obj.Properties = MSProperties_(varargin{:});
            obj.Energy = MSEnergy_(obj.Properties);
            obj.SourceFields = Source_Fields_; 
            obj.State = MSState_;
            obj.EffectiveFields = MSEffectiveFields_(obj);
        end
        
        %% Set Functions
        function obj = Set_Static_Fields(obj,varargin)
            temp = Source_Fields_(varargin{:},'static');
            obj.SourceFields.static = temp.static;
        end
        function obj = Set_Dynamic_Fields(obj,varargin)
            temp = Source_Fields_(varargin{:},'dynamic');
            obj.SourceFields.dynamic = temp.dynamic;
        end
        function obj = Set_State(obj,varargin)
            obj.State = MSState_(obj,varargin{:});
        end
            
        %% update values
        function obj = Update(obj)
            obj.Properties.check_Properties;
            obj.Energy = obj.Energy.Update(obj);
            obj.SourceFields = obj.SourceFields.Update(obj);
            obj.EffectiveFields = obj.EffectiveFields.Update(obj);
        end      
        
    end
    %% Copy function
    methods(Access = protected)
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all four properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the objects
            cpObj.Properties = copy(obj.Properties);
            cpObj.Energy = copy(obj.Energy);
            cpObj.SourceFields = copy(obj.SourceFields);
            cpObj.State = copy(obj.State);
            cpObj.EffectiveFields = copy(obj.EffectiveFields);            
        end
    end
   
    %% Events
    %Add listeners in the future to auto update the particle

    
end


%% Local Functions: Available to the Class only



