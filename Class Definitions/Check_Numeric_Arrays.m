function [ obj ] = Check_Numeric_Arrays(obj,Property_List,err_name)
%CHECK_NUMERIC_ARRAYS Summary of this function goes here
%   Detailed explanation goes here

%ensure all values are numeric and not NaN
sz = zeros(1,numel(Property_List));
for i = 1:numel(Property_List)
    val = obj.(Property_List{i});
    sz(i) = numel(val);
    nanFlag = isnan(val);
    numFlag = ~isnumeric(val);
    if any(nanFlag)
        error('NaN numbers encountered in %s',Property_List{i})
    elseif any(numFlag)
       error('Non-numbers encountered in %s',Property_List{i}) 
    end
end
%make sure all arrays are the same size
switch all(sz == sz(1));
    case 0
        error('Not all %s are the same size',err_name)
    case 1
        
end

end

