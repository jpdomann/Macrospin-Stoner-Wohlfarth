function [particle_num, Xstr, Ystr] = PropertyValue3(input)
%PROPERTYVALUE separate the input cell array into property-value triplets
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