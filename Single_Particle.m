%% Stoner-Wohlfarth Macrospin Particle - Single Particle Example
% John Domann
% 2-1-2017

clear
clc

%Tell matlab where to find the macrospin model
addpath('Class Definitions')

%% Runtime options

freq = 1e9;     %frequency (Hz)
ncycles_dyanamic = 2.0;  %number of cycles
ncycles_static = 1.0;  

nPtsPerCycle_dynamic = 5e3;     %number of time points per cycle
nPtsPerCycle_static = 3e3;
time_dynamic = linspace(0,1/(freq)*ncycles_dyanamic,ncycles_dyanamic*nPtsPerCycle_dynamic);         %times
time_static = linspace(0,1/(freq)*ncycles_static,ncycles_static*nPtsPerCycle_static);         %times

% delT = min(diff(time));
delT = 1e-12;

npts_dynamic = numel(time_dynamic);     %number of time points
npts_static = numel(time_static);     %number of time points

%% Constants 
mu0 = 4*pi*1e-7;

%% Driving fields
%Applied magnetic field (A/m)
H0 = 3e4;
[H1_stat,H2_stat,H3_stat] = MS.Create_Field([H0,0,0],[0 0 0],...
    {'cos','const','const'},nPtsPerCycle_static,npts_static);
[H1_dyn,H2_dyn,H3_dyn] = MS.Create_Field([H0,0,0],[0 0 0],...
    {'cos','const','const'},nPtsPerCycle_dynamic,npts_dynamic);

%Appliedd strain (-)
SMax = 0e-6;
[S1_stat,S2_stat,S3_stat,S4_stat,S5_stat,S6_stat] = MS.Create_Field([SMax,0,0,0,0,0],[0,0,0,0,0,0],...
    {'cos','const','const','const','const','const'},nPtsPerCycle_static,npts_static);
[S1_dyn,S2_dyn,S3_dyn,S4_dyn,S5_dyn,S6_dyn] = MS.Create_Field([SMax,0,0,0,0,0],[0,0,0,0,0,0],...
    {'sin','const','const','const','const','const'},nPtsPerCycle_dynamic,npts_dynamic);

%SOT current (dynamic fields only)
IMax = 0e4;  
[sigma1,sigma2,sigma3] = MS.Create_Field([0,IMax,0],[0 0 0],...
    {'sin','const','const'},nPtsPerCycle_dynamic,npts_dynamic);


%% Create particle
% s = MSParticle_; %with all default properties

% A {'Property','Value'} pair list can be used in MSParticle_ constructor
% Get available properties with:
% prop_list = s.Properties.Property_List; 

%Example to construct a rectangular iron particle
n = MSParticle_('Mat_Name','Nickel','Shape','Ellipse','Dims',[100 80 10]*1e-9,'Location',[0 0 0]);
% n.Properties.Crystal = 'Amorphous';     %and change other properties
% n.Properties.Dims = [100 100 100]*1e-9;   %update the dimensions
% n.Properties.Keb = -mu0*n.Properties.Ms*0.5e4;
n.Properties.Keb = @(s1,s2,s3) -mu0*n.Properties.Ms*0.5e4 * s1;
n.Properties.dir_eb = [1 0 0]; 
n.Properties.Alpha = 5e-1;              %Gilbert Damping
n.Set_State('m',[.1 .1 1]);               %initial magnetization
n.Update;

%% Update Source Fields
n.Set_Static_Fields(...
    'H1',H1_stat,'H2',H2_stat,'H3',H3_stat,...
    'S1',S1_stat,'S2',S2_stat,'S3',S3_stat,...
    'S4',S4_stat,'S5',S5_stat,'S6',S6_stat);
n.Set_Dynamic_Fields(...                    %SOT debug
    'H1',H1_dyn,'H2',H2_dyn,'H3',H3_dyn,...
    'S1',S1_dyn,'S2',S2_dyn,'S3',S3_dyn,...
    'S4',S4_dyn,'S5',S5_dyn,'S6',S6_dyn,...
    'sigma1',sigma1,'sigma2',sigma2,'sigma3',sigma3,...
    't',time_dynamic);

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
%Run MH loops
model.set_static_solver_options('TolX',1e-6,'TolFun',1e-6);
model.set_static_solver_type('segregated');
model.RUN_static;

% model.set_static_solver_type('coupled');
% model.RUN_static;

%% Visualize the results
model.Plots.fig_num = 2;
model.Plots.Particle_Plot(...
    1,'static.H1','Stat_m1',...
    1,'static.H1','Stat_m2',...
    1,'static.H1','Stat_m3');
h1 = model.Plots.plot_handle; %get model handles (figure, axis, lines, and legend)
h1.ax.FontSize = 12;
h1.ax.FontWeight = 'bold';

figure(3)
clf
plot(H1_stat,n.State.Stat_m1)
grid on
[cg,hyst] = MS.centroid(H1_stat(:),n.State.Stat_m1(:));
hold on
plot(cg(1),cg(2),'r*')
text(.1, .65,sprintf('H_{EB} = %2.2e', cg(1)),'Units','Normalized')

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
% model.Plots.make_movie(1,{'dynamic'},'LLG - test movie')

disp('end')


