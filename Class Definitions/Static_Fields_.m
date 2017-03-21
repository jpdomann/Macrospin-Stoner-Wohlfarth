classdef Static_Fields_ < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        H1
        H2
        H3
        S1
        S2
        S3
        S4
        S5
        S6
    end
    properties (Access = public, Constant = true)
       Property_List = {'H1','H2','H3','S1','S2','S3','S4','S5','S6'}; 
    end
    methods
        function obj = Static_Fields_(varargin)
            if nargin <=1 %default constructor
                obj.H1 = 0;
                obj.H2 = 0;
                obj.H3 = 0;
                obj.S1 = 0;
                obj.S2 = 0;
                obj.S3 = 0;
                obj.S4 = 0;
                obj.S5 = 0;
                obj.S6 = 0;
            else
                %separate input into property / value pairs
                [prop, val] = MS.PropertyValue(varargin);
                
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
                    obj.(prop{i}) = val{i};
                end
                
                %Assign properties not described
                NotIncluded = Allowable_Props(~ismember(Allowable_Props,prop));                
                for i = 1:numel(NotIncluded)
                    obj.(NotIncluded{i}) = 0;
                end
                
                %Check MSParticle is in a valid state
                obj = Check_Numeric_Arrays(obj,obj.Property_List,'Static Field Arrays');
            end
            
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
%% Local Functions
