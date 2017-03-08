%% Stoner-Wohlfarth Macrospin Particle - Test Script
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
ncycles = 7;
time = linspace(0,1/(freq)*ncycles,ncycles*nPtsPerCycle);         %times
% delT = min(diff(time));
delT = 1e-12;
npts = numel(time);
%% Driving fields
%Applied magnetic field (A/m)
H0 = 0e6;

% H1 = 0*cos(2*pi*freq*time);
% H2 = 0*sin(2*pi*freq*time);
% H3 = H0*sin(2*pi*freq*time);
 
H1 = linspace(0,0 ,npts);
H2 = linspace(0,0,npts);
H3 = linspace(0,0 ,npts);

% H1 = 0*ones(1,numel(time));
% H2 = 0*ones(1,numel(time));
% H3 = H0*ones(1,numel(time));

%Appliedd strain (-)
SMax = 1000e-6;

S1_1 = zeros(size(H1)); 
S1_2 = SMax*rectpuls(time-3e-9,2e-9); 
S1_3 = SMax*rectpuls(time-4e-9,2e-9); 
S1_4 = SMax*rectpuls(time-5e-9,2e-9); 
S2 = zeros(size(H1)); 
S3 = S2;
S4 = S2; S5 = S2; S6 = S2;

%% Create particle
% s = MSParticle_; %with all default properties

% A 'Property','Value' pair list can be used in MSParticle_ constructor
% prop_list = s.Properties.Property_List;

%Example to construct a rectangular iron particle
n = MSParticle_('Mat_Name','Nickel','Shape','Ellipse','Dims',[1 1 1],'Location',[0 0 0]);
n.Properties.Crystal = 'Amorphous';     %and change other properties
n.Properties.Dims = [100 90 10]*1e-9;   %update the dimensions
n.Properties.Alpha = 1e-1;              %Gilbert Damping
n.Set_State('m',[1 .1 0]);               %initial magnetization
n.Update;

%% Update Source Fields
n.Set_Static_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1_1,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6);
n.Set_Dynamic_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1_1,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6,...
    't',time);

%% Create new particles
%copy to new particle, and change initial magnetization orientation. THe
%source fields applied to 'n' will be copied as well
m=copy(n);
m.Properties.Location = [0 140 0]*1e-9;
m.Set_Dynamic_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1_2,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6,...
    't',time);
m.Update;

p=copy(n);
p.Set_State('m',[-1 .1 0]); 
p.Set_Dynamic_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1_3,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6,...
    't',time);p.Properties.Location = [0 2*140 0]*1e-9;
p.Update;

q=copy(m);
q.Set_Dynamic_Fields(...
    'H1',H1,'H2',H2,'H3',H3,...
    'S1',S1_4,'S2',S2,'S3',S3,...
    'S4',S4,'S5',S5,'S6',S6,...
    't',time);
q.Properties.Location = [0 3*140 0]*1e-9;
q.Update;

%% Create MS Model
particles = {n,m,p,q};
% particles = {n,m};
% particles = {n};
nparts = numel(particles);
model = MSModel_(particles);
%update ode options 
% detT = 1e-11;
model.set_ode_options('MaxStep',delT,'InitialStep',delT*1e-3');   %solvers: 'ode45','ode23','ode113','ode15s'
model.r_cutoff = inf; %0 - no dipole coupling, inf - couples all particles

%% Check initial setup
model.Plots.fig_num = 1;
model.Plots.Plot_Initial_State(1:4);

%% Dynamic MH Loops
model.set_ode_options('solver','ode15s');
model.r_cutoff = inf; %dipole couple all particles
model.RUN_LLG;

%Sample times:
%ode45: 55sec   %ode23: 33sec    %ode113: 23sec    %ode15s: 15sec 

%% Visualize the results
% For possible values, see:
% model.Plots.Source_field_names
% model.Plots.State_names

model.Plots.fig_num = 3;
model.Plots.Particle_Plot(...
    1,'dynamic.t','Dyn_m1',...
    1,'dynamic.t','Dyn_m2',...
    2,'dynamic.t','Dyn_m1',...
    2,'dynamic.t','Dyn_m2',...
    3,'dynamic.t','Dyn_m1',...
    3,'dynamic.t','Dyn_m2',...
    4,'dynamic.t','Dyn_m1',...
    4,'dynamic.t','Dyn_m2');
h2 = model.Plots.plot_handle;
for i = 1:numel(h2.lines); h2.lines(i).LineWidth = 2; end
h2.ax.FontWeight = 'bold';

%% Animate
model.Plots.fig_num = 1;
model.Plots.Animate_Spins([1,2,3,4],{'dynamic'})

%% Save Movie
model.Plots.make_movie([1,2,3,4],{'dynamic'},'bennet-clocking-2d')

disp('end')


%% Old code
% %% Static 
% % figure(1)
% % clf
% % start=1;
% % markers = {'s','o','<','>'};
% % for i = 1:nparts
% %     theta = atan2(sqrt((model.particles{i}.State.Stat_m1).^2 +(model.particles{i}.State.Stat_m2).^2),model.particles{i}.State.Stat_m3);
% %     phi = atan2(model.particles{i}.State.Stat_m2,model.particles{i}.State.Stat_m1);
% %     
% %     ht = plot(H1(start:end), theta(start:end)*180/pi,'-');
% %     hold all
% %     hp = plot(H1(start:end),phi(start:end)*180/pi,'-');
% %     
% %     ht.LineWidth = 5-i;
% %     ht.MarkerSize = 5-i;
% %     hp.LineWidth = 5-i;
% %     
% % end
% % grid on
% 
% figure(3)
% clf
% Hp = H3;
% for i = 1:nparts
%     
%     m1 = model.particles{i}.State.Stat_m1;
%     m2 = model.particles{i}.State.Stat_m2;
%     m3 = model.particles{i}.State.Stat_m3;
%     %     comet3(m1,m2,m3);
%     comet(Hp,m3)
%     h = plot(Hp,m1,Hp,m2,Hp,m3);    
%     legend('m1','m2','m3')
% %     h.LineWidth = 5-i;
% %     h.MarkerSize = 5-i;   
% 
% %     axis([-1 1 -1 1 -1 1])
%     grid on
% 
%     drawnow      
%     pause(1)
% end
% 
% %% Dynamic
% % figure(4)
% % clf
% % start=1;
% % for i = 1:nparts
% %     theta = atan2(sqrt((model.particles{i}.State.Dyn_m1).^2 +(model.particles{i}.State.Dyn_m2).^2),model.particles{i}.State.Dyn_m3);
% %     phi = atan2(model.particles{i}.State.Dyn_m2,model.particles{i}.State.Dyn_m1);
% %     
% %     ht = plot(H1(start:end), theta(start:end)*180/pi,'-');
% %     hold all
% %     hp = plot(H1(start:end),phi(start:end)*180/pi,['-',markers{i}]);
% %     
% %     ht.LineWidth = 5-i;
% %     ht.MarkerSize = 5-i;
% % end
% % grid on
% 
% 
% figure(6)
% clf
% Hp = H2;
% for i = 1:nparts
%     
%     m1 = model.particles{i}.State.Dyn_m1;
%     m2 = model.particles{i}.State.Dyn_m2;
%     m3 = model.particles{i}.State.Dyn_m3;
%     %     comet3(m1,m2,m3);
%     comet(Hp,m3)
%     h = plot(Hp,m1,Hp,m2,Hp,m3);    
%     legend('m1','m2','m3')
% %     h.LineWidth = 5-i;
% %     h.MarkerSize = 5-i;   
% 
% %     axis([-1 1 -1 1 -1 1])
%     grid on
% 
%     drawnow      
%     pause(1)
% end
% 
% 
% for i = 1:nparts
%     figure(6+i)
%     if i == 1; clf; end    
%     m1 = model.particles{i}.State.Dyn_m1;
%     m2 = model.particles{i}.State.Dyn_m2;
%     m3 = model.particles{i}.State.Dyn_m3;
% %     comet3(m1,m2,m3,1e-2);
%     plot3(m1,m2,m3,'-',m1(1),m2(1),m3(1),'*',m1(end),m2(end),m3(end),'>')
%     hold all
%     axis([-1 1 -1 1 -1 1])
%     grid on
%     drawnow      
% end
% 
% 
% figure (8)
% clf
% for i = 1:nparts
%     m1 = model.particles{i}.State.Dyn_m1;
%     m2 = model.particles{i}.State.Dyn_m2;
%     m3 = model.particles{i}.State.Dyn_m3;
%     
%     plot(time,[m1,m2,m3])
%     hold all
%     
% end
% legend('m1','m2','m3','n1','n2','n3')
% 
% %% Dynamic II
% figure(9)
% clf
% subplot(2,1,2)
% h1 = plot3(0,0,0,'-');
% hold on
% h2 = plot3(0,0,0,'*r');
% axis([-1 1 -1 1 -1 1])
% grid on
% 
% m1 = model.particles{2}.State.Dyn_m1;
% m2 = model.particles{2}.State.Dyn_m2;
% m3 = model.particles{2}.State.Dyn_m3;
% 
% for i = 1:5:numel(m1)
%     h1.XData = m1(1:i);
%     h1.YData = m2(1:i);
%     h1.ZData = m3(1:i);
%     
%     h2.XData = m1(i);
%     h2.YData = m2(i);
%     h2.ZData = m3(i);
%     
%     drawnow
% end
% 
% 
%    