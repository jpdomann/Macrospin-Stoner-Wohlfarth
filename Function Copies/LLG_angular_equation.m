function  [dangdt,m,Heff] = LLG_angular_equation(t,ang_in,model,particles,nparts,time,pb)
%LLG_EQUATION Summary of this function goes here
%   Input:
%       t - time
%       ang_in - magnetization orientation (theta,phi)
%       model - handle to model
%       particles - handle to particles in model
%       nparts - number of particles
%       time - entire time array
%       pb - progress bar
%   Output:
%       dangdt - time rate of change of magnetization orientation angles
%       m - updated magnetization orientation
%       Heff - effective magnetic field

%% Format input vectors
%resize
ang = reshape(ang_in,2,nparts);

%compute magnetization components (spin)
m = zeros(3,nparts);
m(1,:) = sin(ang(1,:)).*cos(ang(2,:));
m(2,:) = sin(ang(1,:)).*sin(ang(2,:));
m(3,:) = cos(ang(1,:));

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
%store factors for each particle
F1 = zeros(3,nparts);
F2 = F1;
for i = 1:nparts
    F1(1:3,i) = particles{i}.Properties.Factor_1 ;
    F2(1:3,i) = particles{i}.Properties.Factor_2 ;
end

%Vector form
dmdt = - F1 .* cross(m,Heff) + F2 .* cross(m,cross(m,Heff)) ;

%% Convert to dang/dt
%convert m_dot to theta_dot and phi_dot
dangdt = zeros(2,nparts);
%theta_dot
dangdt(1,:) = sec(ang(1,:)).*cos(ang(2,:)).*dmdt(1,:) + ...
    sec(ang(1,:)).*sin(ang(2,:)).*dmdt(2,:);
%phi_dot
dangdt(2,:) = -csc(ang(1,:)).*sin(ang(2,:)).*dmdt(1,:) + ...
    csc(ang(1,:)).*cos(ang(2,:)).*dmdt(2,:);


%% update model
model.time{model.count} = t;
model.m{model.count} = m;
model.ang{model.count} = ang;
model.Heff{model.count} = Heff;
model.count = model.count + 1;

%% output variables
dangdt = dangdt(:);
m = m(:);
Heff = Heff(:);

%% Update progress bar
progress = floor((t-time(1)) / (time(end) - time(1))*1e3);
diff = progress - model.progress;
switch diff > 0
    case 1
        for i = 1:diff
            pb.progress;
            model.progress = model.progress + 1;
        end
    case 0        
end

end

