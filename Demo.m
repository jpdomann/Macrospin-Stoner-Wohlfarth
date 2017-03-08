%% Stoner-Wohlfarth Macrospin Particle - Demo Script
% John Domann
% 3-7-2017
clear
clc

%Tell matlab where to find the macrospin model
addpath('Class Definitions')

%% Setup time points
%time arrays
nPtsPerCycle_dynamic = 1e3;     	%dynamic - number of time points
nPtsPerCycle_static = 2e2;     	%static - number of time points
freq = 1e9;
ncycles = 2;
time_dynamic = linspace(0,1/(freq)*ncycles,ncycles*nPtsPerCycle_dynamic);         %time
time_static = linspace(0,1/(freq)*ncycles,ncycles*nPtsPerCycle_static);         %time

delT = min(min(diff(time_dynamic)), 1e-12); %minimum time stepping in LLG ODE Solver
npts_dynamic = numel(time_dynamic);
npts_static = numel(time_static);

%% Driving fields
%Applied magnetic field (A/m)
H0 = 1e6;

%Fields
H1 = @(time) 0*cos(2*pi*freq*time);             %sinusoidal profile
H2 = @(time) 0*sin(2*pi*freq*time);
H3 = @(time) H0*cos(2*pi*freq*time);
 
% H1 = @(time) linspace(0,0 ,numel(time));      %linear ramp
% H2 = @(time) linspace(0,0,numel(time));
% H3 = @(time) linspace(H0,0 ,numel(time));

% H1 = @(time) 0*ones(1,numel(time));           %contant fields
% H2 = @(time) 0*ones(1,numel(time));
% H3 = @(time) H0*ones(1,numel(time));

%Appliedd strain (-)
SMax = 0e-6;

S1 = @(time) SMax*ones(size(H1(time))); 
S2 = @(time) S1(time); 
S3 = @(time) S2(time);
S4 = @(time) S2(time); S5 = @(time) S2(time); S6 = @(time) S2(time);

%% Create particle
% s = MSParticle_; %with all default properties

% A 'Property','Value' pair list can be used in MSParticle_ constructor
% prop_list = s.Properties.Property_List;

%Example to construct a rectangular iron particle
n = MSParticle_('Mat_Name','Nickel','Shape','Rectangle','Dims',[5 1 2],'Location',[0 0 0]);

%Values can be updated after creation
n.Properties.Shape='Ellipse';       %you can update shape after creation
n.Properties.Crystal = 'Amorphous'; %and change other properties
n.Properties.Dims = [100 100 100]*1e-9;        %or update the dimensions
n.Properties.Alpha = 1e-1;

%Set_State will normalize all inputs to a magnitude of 1, so inputting
%crystal directions like [1 1 1] works just fine
n.Set_State('m1',0,'m2',0.1,'m3',1);    % initial magnetic orientation
n.Set_State('m2',0.1,'m3',1);           % you only specify non-zero components
n.Set_State('m',[1 1 .1]);             % or specify all components with a vector

%But make sure to update the particle! (TRY REPLACING WITH LISTENER PROP.)
n.Update;

%% Update Source Fields
n.Set_Static_Fields(...
    'H1',H1(time_static),'H2',H2(time_static),'H3',H3(time_static),...
    'S1',S1(time_static),'S2',S2(time_static),'S3',S3(time_static),...
    'S4',S4(time_static),'S5',S5(time_static),'S6',S6(time_static));

n.Set_Dynamic_Fields(...
    'H1',H1(time_dynamic),'H2',H2(time_dynamic),'H3',H3(time_dynamic),...
    'S1',S1(time_dynamic),'S2',S2(time_dynamic),'S3',S3(time_dynamic),...
    'S4',S4(time_dynamic),'S5',S5(time_dynamic),'S6',S6(time_dynamic),...
    't',time_dynamic);

%% Create new particles
%copy to new particle, and change initial magnetization orientation. THe
%source fields applied to 'n' will be copied as well
m=copy(n);
m.Set_State('m1',-1,'m3',.1); 
m.Properties.Location = [0 140 0]*1e-9;
% m.Properties.Alpha = 1e-1;
m.Update;

%% Create MS Model
particles = {n,m};              %group particles together
nparts = numel(particles);      %store number of particles
model = MSModel_(particles);    %create a macrospin model with the particles

%update ode options 
model.set_ode_options('MaxStep',delT,'InitialStep',delT*1e-3');   %solvers: 'ode45','ode23','ode113','ode15s'
model.r_cutoff = inf;            %dipole coupling radius: 0 = no dipole coupling, inf = couples all particles

%% Check Initial State
%model.Plots.Plot_Initial_State([Particle Numbers]);
% model.Plots.Plot_Initial_State([1]);          %shows only particle 1
% model.Plots.Plot_Initial_State([2]);          %shows only particle 2
model.Plots.Plot_Initial_State([1:nparts]);     %show all particles

%% Static MH Loops
%Run MH loops
model.set_static_solver_options('TolX',1e-6,'TolFun',1e-6)
% model.set_static_solver_type('segregated'); %STILL NEEDS WORK
model.set_static_solver_type('coupled');


model.RUN_static;

%% Visualize the static results
% For possible values, see:
% model.Plots.Source_field_names
% model.Plots.State_names
model.Plots.fig_num = 1;
model.Plots.Particle_Plot(...   
    1,'static.H3','Stat_m1',...    
    1,'static.H3','Stat_m2',...    
    1,'static.H3','Stat_m3',...    
    2,'static.H3','Stat_m1',...
    2,'static.H3','Stat_m2',...
    2,'static.H3','Stat_m3');
h1 = model.Plots.plot_handle; %get model handles (figure, axis, lines, and legend)
h1.ax(1).LineWidth = 1.5;         %so you can refine how the plot looks
h1.ax(1).FontSize = 12;
h1.ax(1).FontWeight = 'bold';
for i = 1:numel(h1.lines); h1.lines(i).LineWidth = 2; end

%% Animate Static
model.Plots.fig_num = 2;
model.Plots.Animate_Spins([1:nparts],{'static'})

%% Dynamic MH Loops
model.set_ode_options('solver','ode15s');
model.r_cutoff = inf; %dipole couple all particles
model.RUN_LLG;

%Sample times:
%ode45: 55sec 
%ode23: 33sec
%ode113: 23sec
%ode15s: 15sec    ----      RECOMENDED SOLVER

%% Visualize the dynamic results
% For possible values, see:
% model.Plots.Source_field_names
% model.Plots.State_names
model.Plots.fig_num = 3;
model.Plots.Particle_Plot(...
    1,'dynamic.t','Dyn_m1',...%
    1,'dynamic.t','Dyn_m2',...
    1,'dynamic.t','Dyn_m3',...
    2,'dynamic.t','Dyn_m1',...%
    2,'dynamic.t','Dyn_m2',...
    2,'dynamic.t','Dyn_m3');
h2 = model.Plots.plot_handle;
for i = 1:numel(h2.lines); h2.lines(i).LineWidth = 2; end

%% Animate Dynamic
model.Plots.fig_num = 4;
model.Plots.Animate_Spins([1:nparts],{'dynamic'})

%% Save movies
model.Plots.fig_num = 2;
model.Plots.make_movie([1:nparts],{'static'},'Movies\demo-movie-static')
model.Plots.fig_num = 4;
model.Plots.make_movie([1:nparts],{'dynamic'},'Movies\demo-movie-dynamic')


disp('end')

