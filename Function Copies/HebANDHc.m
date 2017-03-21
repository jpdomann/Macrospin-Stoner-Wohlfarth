function [ Hc,Heb, Mc, Meb ] = HebANDHc( H,M )
%HEBANDHC Summary of this function goes here
%   Detailed explanation goes here

%Break into rising and falling segments
I_min = find(H == min(H));
I_min = I_min(1);

if I_min == 1
    I_min = find(H == max(H));
    I_min = I_min(1);
    Hrise = H(1:I_min);
    Mrise = M(1:I_min);
    Hfall = H(I_min:end);
    Mfall = M(I_min:end);
else
    Hfall = H(1:I_min);
    Mfall = M(1:I_min);
    Hrise = H(I_min:end);
    Mrise = M(I_min:end);
end

%Ensure monotonic behavior of data (add / subtract eps 
Hfall = makeMonotonic(Hfall);
Mfall = makeMonotonic(Mfall);
Hrise = makeMonotonic(Hrise);
Mrise = makeMonotonic(Mrise);

%Fit an interpolation curve for each region
fit = 'linear';
n=1e4;
hfall = linspace(Hfall(1),Hfall(end),n);
mfall = interp1(Hfall,Mfall,hfall,fit,'extrap');

hrise = linspace(Hrise(1),Hrise(end),n);
mrise = interp1(Hrise,Mrise,hrise,fit,'extrap');

%Find crossing locations
Ifall = find(mfall<0,1);
Irise = find(mrise>0,1);

%Coercive field
Hc(1) = hfall(Ifall);
Hc(2) = hrise(Irise);
Mc(1) = mfall(Ifall);
Mc(2) = mrise(Irise);
%Exchange Bias
Heb = mean(Hc);
Meb = mean(Mc);

% %plot
% plot(Hfall,Mfall,'r',Hrise,Mrise,'b')
% hold on
% plot(hfall,mfall,'b.',hrise,mrise,'r.')
% plot(Hc,Mc,'g*',Heb,Meb,'g*')
% grid on

end

