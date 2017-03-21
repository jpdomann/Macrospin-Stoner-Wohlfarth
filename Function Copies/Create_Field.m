function [ varargout ] = Create_Field( amplitude, ave, shape, pts_per_cycle, npts )
%CREATE_FIELD has a set of built in waveforms that can be used to create
%the various input fields. 
%   Inputs: 
%       ave - average value of waveform
%       amplitude - amplitude of waveformf (half peak to peak)
%       shape - sin, cos, ramp, or const
%       pts_per_cycle - number of data points per cycle
%       npts - overall length of field ( note: the number of cycles is
%       npts/pts_per_cycle)
%   Output:
%       n arrays of size (1 x npts). n is the length of the input amplitdue
%       array
%   Notes: If ave or shape don't have the same number of elements as
%   amplitdue, then they must have only one value, which will be used for
%   all arrays.

if nargout ~= numel(amplitude)
    error('Number of outputs must be the same as the number of amplitudes');
end

%if only one shape is entered, apply for all fields
switch numel(amplitude) ~= numel(shape) || numel(amplitude) ~= numel(ave)
    case 1
        switch numel(shape)
            case 1
                shape = repmat(shape,numel(amplitude),1);
            otherwise
                error('Number of input shapes must equal number of input amplitudes, or be equal to 1')
        end
        
        switch numel(ave)
            case 1
                ave = repmat(ave,numel(amplitude),1);
            otherwise
                error('Number of input ave must equal number of input amplitudes, or be equal to 1')
        end
end
varargout = cell(1, nargout);

for i = 1:numel(amplitude)      
    A = amplitude(i);   %amplitude of given wave
    switch shape{i} %create fields
        case 'sin'
            varargout{i} = A*sin([1:npts]*2*pi/pts_per_cycle) + ave(i);
        case 'cos'
            varargout{i} = A*cos([1:npts]*2*pi/pts_per_cycle) + ave(i);
        case 'ramp'
            temp = A*(linspace(0,1,pts_per_cycle) - 1/2) + ave(i);
            for j = 1:ceil(npts/pts_per_cycle)
                temp = [temp(:) temp(:)];
            end
            varargout{i} = temp(1:npts);
        case 'const'
            varargout{i} = A*ones(1,npts);
    end       
end

end