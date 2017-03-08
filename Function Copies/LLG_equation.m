function  [dmdt,m,Heff] = LLG_equation(t,m_in,model,particles,nparts,time)
%LLG_EQUATION Summary of this function goes here
%   Input:
%       t - time
%       m_in - magnetization orientation (spin)
%       model - handle to model
%       particles - handle to particles in model
%       nparts - number of particles
%       time - entire time array
%   Output:
%       dmdt - time rate of change of magnetization orientation
%       m - updated magnetization orientation
%       Heff - effective magnetic field

%% Format input vectors
%resize 
m = reshape(m_in,3,nparts);

% Normalize input state (abs(m)=1)
normM = sum(m.^2);
m = bsxfun(@rdivide,m,normM);

%% Get effective fields
%linearly interpolate between given points 
Heff = zeros(3,nparts);
for i = 1:nparts
    f = particles{i}.SourceFields.dynamic;
    fint = @(name,t) interp1(time,f.(name),t);
    Heff(1,i) = particles{i}.EffectiveFields.H1_total(...
        fint('H1',t),fint('H2',t),fint('H3',t),...
        m(1,i),m(2,i),m(3,i),...
        fint('S1',t),fint('S2',t),fint('S3',t),fint('S4',t),fint('S5',t),fint('S6',t));
    Heff(2,i) = particles{i}.EffectiveFields.H2_total(...
        fint('H1',t),fint('H2',t),fint('H3',t),...
        m(1,i),m(2,i),m(3,i),...
        fint('S1',t),fint('S2',t),fint('S3',t),fint('S4',t),fint('S5',t),fint('S6',t));
    Heff(3,i) = particles{i}.EffectiveFields.H3_total(...
        fint('H1',t),fint('H2',t),fint('H3',t),...
        m(1,i),m(2,i),m(3,i),...
        fint('S1',t),fint('S2',t),fint('S3',t),fint('S4',t),fint('S5',t),fint('S6',t));
    
end

%% Compute dm/dt
%Vector form
dmdt = - particles{i}.Properties.Factor_1 .* cross(m,Heff)...
    + particles{i}.Properties.Factor_1 .* cross(m,cross(m,Heff)) ;

%% update model

end

