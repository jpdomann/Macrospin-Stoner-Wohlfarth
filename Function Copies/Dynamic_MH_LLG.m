function [ model ] = Dynamic_MH_LLG( model )
%DYNAMIC_MH_LLG Takes a SWModel_ and a list of the 
%   The f

%% Create anonymous ODE function to be solved

%model particles
particles = model.particles;    %particles stored in model
nparts = numel(particles);      %number of particles
% npts = numel(particles{1}.SourceFields.dynamic.H1); %number of field points

%time interval
time = particles{1}.SourceFields.dynamic.t; %specify specfic time points
% tspan = [min(time), max(time)]; %use to specify endpoint only

%initial magnetization state (3 x n) and angle state (2 x n) 
m0 = zeros(3,nparts); %(m1,m2,m3)
ang0 = zeros(2,nparts);%(theta, phi)
for i = 1:nparts
    m0(1,i) = particles{i}.State.m1;                            %m1
    m0(2,i) = particles{i}.State.m2;                            %m2
    m0(3,i) = particles{i}.State.m3;                            %m3
    ang0(1,i) = atan2(sqrt(m0(1,i)^2 + m0(2,i)^2) , m0(3,i) );  %theta
    ang0(2,i) = atan2(m0(2,i), m0(1,i) );                       %phi
end
% model.m{model.count} = m0;

%progress bar
p = ProgressBar(1e3);

% anonymous ODE to be solved
dang_dt = @(t,ang) LLG_angular_equation(t,ang,model,particles, nparts,time,p);

%% Solve LLG Equation
options = odeset('RelTol',1e-6,'AbsTol',1e-9,'NormControl','on','InitialStep',1e-14,'MaxStep',1e-10);

[t_solved,ang_solved]  = ode45(dang_dt, time, ang0(:),options);
p.stop;

%ensure t_solved = time
 switch all(t_solved(:) == time(:));
    case 0
        warning('ode solver output different times than those input')   
end

%% Update Particles with solution
for i = 0:nparts-1    
    particles{i+1}.State.set_state(sin(ang_solved(:,2*i+1)).*cos(ang_solved(:,2*i+2)),'Dyn_m1');
    particles{i+1}.State.set_state(sin(ang_solved(:,2*i+1)).*sin(ang_solved(:,2*i+2)),'Dyn_m2');
    particles{i+1}.State.set_state(cos(ang_solved(:,2*i+1)),'Dyn_m3');    
end


end

