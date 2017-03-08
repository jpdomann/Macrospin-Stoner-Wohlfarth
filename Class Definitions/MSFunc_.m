classdef MSFunc_ < handle
    %MSFUNC_ Is a utility class that holds useful functions to be accessed
    %by multiple other classes
    
    properties
        prop1
    end
    
    methods
        function obj = MSFunc_(obj)
            
        end
        
        %% PropertyValue         
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
        
        
    end
    
end

