function [ varargout ] = Create_Field( amplitude, shape, pts_per_cycle, npts )
%CREATE_FIELD Summary of this function goes here
%   Detailed explanation goes here

if nargout ~= numel(amplitude)
    error('Number of outputs must be the same as the number of amplitudes');
end

%if only one shape is entered, apply for all fields
switch numel(amplitude) ~= numel(shape)
    case 1
        switch numel(shape)
            case 1
                shape = repmat(shape,numel(amplitude),1);
            otherwise
                error('Number of input shapes must equal number of input amplitudes, or be equal to 1')
        end
end
varargout = cell(1, nargout);

for i = 1:numel(amplitude)
    A = amplitude(i);
    switch shape{i}
        case 'sin'
            varargout{i} = A*sin([1:npts]*2*pi/pts_per_cycle);
        case 'cos'
            varargout{i} = A*cos([1:npts]*2*pi/pts_per_cycle);
        case 'ramp'
            temp = A*linspace(0,1,pts_per_cycle);
            for j = 1:ceil(npts/pts_per_cycle)
                temp = [temp(:) temp(:)];
            end
            varargout{i} = temp(1:npts);
        case 'const'
            varargout{i} = A*ones(1,npts);
    end
    
end

end
% H1 = H0*sin(2*pi*freq*time_dynamic);   %sinusoidal
% H2 = 0*sin(2*pi*freq*time_dynamic);
% H3 = 0*cos(2*pi*freq*time_dynamic);

% H1 = linspace(0,H0 ,npts);       %ramp
% H2 = linspace(0,0,npts);
% H3 = linspace(0,0 ,npts);

% H1 = 0*ones(1,numel(time));   %constant
% H2 = 0*ones(1,numel(time));
% H3 = H0*ones(1,numel(time));