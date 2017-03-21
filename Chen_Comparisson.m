%% Stoner-Wohlfarth Macrospin Particle - Chen Comparisson
% John Domann
% 3-21-2017
% This script replicates the results of Chet et. al. "Angular dependence of
% exchange bias and magnetization reversal controlled by electric field
% induced competing anisotropies"

clear
clc

%Tell matlab where to find the macrospin model
addpath('Class Definitions')

%% Runtime options

freq = 1e9;     %frequency (Hz)
ncycles_static = 1.0;   %number of cycles  
nPtsPerCycle_static = 1e3;          %number of time points per cycle
time = linspace(0,1/(freq)*ncycles_static,ncycles_static*nPtsPerCycle_static);	%time
delT = min(diff(time));
npts_static = numel(time);     %number of time points

%% Constants
mu0 = 4*pi*1e-7;
AM2Oe = 1/79.577; %multiply A/m by AM2Oe to get Oe

%angles to test
theta = linspace(0, pi/2,50);

results.particles = cell(1,numel(theta));
results.Heb1 = zeros(1,numel(theta));
results.Heb1_Oe = zeros(1,numel(theta));
results.Heb2 = zeros(1,numel(theta));
results.Heb2_Oe = zeros(1,numel(theta));

for i = 1:numel(theta)
    
    %% Driving fields
    %Applied magnetic field (A/m)
    H0 = 1.5915e4*3/2;  %200 Oe
    [H1_stat,H2_stat,H3_stat] = MS.Create_Field(H0*[cos(theta(i)),sin(theta(i)),0],[0 0 0],...
        {'cos','cos','const'},nPtsPerCycle_static,npts_static);
    
    %Appliedd strain (-)
    %strain will be introduced with an effective uniaxial anisotropy energy
    %instead of the traditional cubic anisotropy expression (to keep it the
    %same as the chen paper)
    [S1_stat,S2_stat,S3_stat,S4_stat,S5_stat,S6_stat] = MS.Create_Field([0,0,0,0,0,0],[0,0,0,0,0,0],...
        {'const','const','const','const','const','const'},nPtsPerCycle_static,npts_static);
    
    %% Create particle
    % construct particle
    n = MSParticle_('Mat_Name','CoFeB','Shape','Ellipse','Dims',[1000 1000 55]*1e-9,'Location',[0 0 0]);
    %update properties based on chen paper
    n.Properties.Crystal = 'Amorphous';     %no crystalline anisotropy
    n.Properties.Ms = 1.11e6;               %update Ms 
    n.Properties.Ks = 0;                    %turn off PMA
    n.Properties.B_me = [0 0];                  %turn of normal magnetoelasticity
    
    %exchange bias
    n.Properties.Keb = -mu0*n.Properties.Ms*2796;   %chen = 3900 J/m^3
    % n.Properties.Keb = @(s1,s2,s3) -mu0*n.Properties.Ms*0.5e4 .* s1;
    n.Properties.dir_eb = [1 0 0];   
    
    %uniaxial anisotropy ( static and stress induced)
    S1Max = -846e-6;
    S2Max = 210e-6;
    Ku1 = 2900;     %J/m^3 - due to bias field (from chen)
    Ku2 = -3/2*(30e-6)*(160e9)*(S2Max - S1Max);     %due to stress ~7600 J/m^3
    n.Properties.Kuni = Ku1;                %just static anisotropy
%     n.Properties.Kuni = Ku1 + Ku2;            %static and stress
    n.Properties.dir_uni = [1 0 0]; 
    
    n.Properties.Alpha = 5e-1;              %Gilbert Damping
    n.Set_State('m',[.1 .1 1]);               %initial magnetization
    n.Update;
    
    %% Update Source Fields
    n.Set_Static_Fields(...
        'H1',H1_stat,'H2',H2_stat,'H3',H3_stat,...
        'S1',S1_stat,'S2',S2_stat,'S3',S3_stat,...
        'S4',S4_stat,'S5',S5_stat,'S6',S6_stat);
 
    %% Create MS Model
    particles = {n};
    nparts = numel(particles);
    model = MSModel_(particles);
    
    %% Static MH Loops
    %Run MH loops
    model.set_static_solver_options('TolX',1e-8,'TolFun',1e-8);
    %     model.set_static_solver_type('segregated');
    %     model.RUN_static;
    tic
    model.set_static_solver_type('coupled');
    model.RUN_static;
    t(i) = toc;
    
    %% Visualize the results
    %projection in measurement direction
    H_theta = H1_stat.*cos(theta(i)) + H2_stat.*sin(theta(i));
    m_theta = n.State.Stat_m1.*cos(theta(i)) + n.State.Stat_m2.*sin(theta(i));
          
    figure(1)
    clf
    plot(H_theta*AM2Oe,m_theta)
    grid on
    [cg,hyst] = MS.centroid(H_theta(:),m_theta(:));
    [Hc,Heb, Mc, Meb] = MS.HebANDHc(H_theta(:),m_theta(:));
    
    hold on
%     plot(cg(1)*AM2Oe,cg(2),'r*')
    plot(Heb*AM2Oe,Meb,'go')
    text(.1, .65,sprintf('H_{EB} = %2.2f Oe', cg(1)*AM2Oe ),'Units','Normalized')
    text(.1, .55,sprintf('H_{EB} = %2.2f Oe', Heb*AM2Oe ),'Units','Normalized')
    xlabel('H (Oe)')
    ylabel('M/M_s (-)')
        
    %% Store results
    results.Heb1(i) = cg(1);
    results.Heb2(i) = Heb;
    results.Heb1_Oe(i) = cg(1)*AM2Oe;
    results.Heb2_Oe(i) = Heb*AM2Oe;
    results.particles{i} = copy(n);
      
    figure(2)
    clf
    %plot(theta*180/pi, -results.Heb1_Oe,'r*-',theta*180/pi, -results.Heb2_Oe,'g*-')
    plot(theta*180/pi, -results.Heb2_Oe,'g*-')
    xlabel('\Theta (degrees)')
    ylabel('H_{eb} (Oe)')
    grid on
    drawnow
    
    clear n
end
results.theta = theta;

%% Save Results
save('results-no-strain','results')


