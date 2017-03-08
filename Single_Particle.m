%% Stoner-Wohlfarth Macrospin Particle - Single Particle Example
% John Domann
% 2-1-2017

clear
clc

%Tell matlab where to find the macrospin model
addpath('Class Definitions')

%% Runtime options
%time array
nPtsPerCycle = 5e3;                           %number of time points
freq = 1e9;
ncycles = 1.0;
time = linspace(0,1/(freq)*ncycles,ncycles*nPtsPerCycle);         %times
% delT = min(diff(time));
delT = 1e-12;
npts = numel(time);

%% Driving fields
%Applied magnetic field (A/m)
H0 = 1e6;

% H1 = 0*cos(2*pi*freq*time);   %sinusoidal
% H2 = 0*sin(2*pi*freq*time);
% H3 = H0*cos(2*pi*freq*time);
 
H1 = linspace(0,0 ,npts);       %ramp
H2 = linspace(0,0,npts);
H3 = linspace(0,H0 ,npts);

% H1 = 0*ones(1,numel(time));   %constant
% H2 = 0*ones(1,numel(time));
% H3 = H0*ones(1,numel(time));

%Appliedd strain (-)
SMax = 0e-6;

S1 = SMax*rectpuls(time-2e-9,2e-9);
S2 = zeros(size(H1)); 
S3 = S2;
S4 = S2; S5 = S2; S6 = S2;

%% Create particle
% s = MSParticle_; %with all default properties

% A {'Property','Value'} pair list can be used in MSParticle_ constructor
% Get available properties with:
% prop_list = s.Properties.Property_List; 

%Example to construct a rectangular iron particle
n = MSParticle_('Mat_Name','Nickel','Shape','Ellipse','Dims',[1 1 1],'Location',[0 0 0]);
n.Properties.Crystal = 'Amorphous';     %and change other properties
n.Properties.Dims = [100 100 100]*1e-9;   %update the dimensions
n.Properties.Alpha = 1e-2;              %Gilbert Damping
n.Set_State('m',[-1 1 1]);               %initial magnetization
n.Update;

%% Update Source Fields
n.Set_Static_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6);
n.Set_Dynamic_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6,...
    't',time);

%% Create MS Model
particles = {n};
nparts = numel(particles);
model = MSModel_(particles);
%update ode options 
% detT = 1e-11;
model.set_ode_options('MaxStep',delT,'InitialStep',delT*1e-3);   %solvers: 'ode45','ode23','ode113','ode15s'
model.r_cutoff = inf; %0 - no dipole coupling, inf - couples all particles

%% Check initial setup
model.Plots.fig_num = 2;
model.Plots.Plot_Initial_State(1);

%% Static MH Loops
% %Run MH loops
% model.set_static_solver_options('TolX',1e-6,'TolFun',1e-6)
% model.set_static_solver_type('segregated');
% model.RUN_static;

% model.set_static_solver_type('coupled');
% model.RUN_static;

%% Dynamic MH Loops
model.set_ode_options('solver','ode15s');
model.r_cutoff = inf; %dipole couple all particles
model.RUN_LLG;

%Sample times:
%ode45: 55sec   %ode23: 33sec    %ode113: 23sec    %ode15s: 15sec 

%% Visualize the results
model.Plots.fig_num = 2;
model.Plots.Particle_Plot(...
    1,'dynamic.H3','Dyn_m1',...
    1,'dynamic.H3','Dyn_m2',...
    1,'dynamic.H3','Dyn_m3');
h1 = model.Plots.plot_handle; %get model handles (figure, axis, lines, and legend)
h1.ax.FontSize = 12;
h1.ax.FontWeight = 'bold';

model.Plots.fig_num = 3;
model.Plots.Particle_Plot(...
    1,'dynamic.t','Dyn_m1',...%1,'dynamic.t','Dyn_m2',...
    1,'dynamic.t','Dyn_m2',...
    1,'dynamic.t','Dyn_m3');
h2 = model.Plots.plot_handle;
for i = 1:numel(h2.lines); h2.lines(i).LineWidth = 2; end
h2.ax.FontWeight = 'bold';

%% Animate
model.Plots.fig_num = 1;
model.Plots.Animate_Spins(1,{'dynamic'})

%% Save Movie
model.Plots.make_movie(1,{'dynamic'},'LLG - test movie')

disp('end')


