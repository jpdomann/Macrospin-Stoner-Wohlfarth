classdef MSModel_<handle
    %MSPARTICLE Assembles collection of MSParticles_ into a single model.
    %It performs basic checks (like making sure all applied field arrays
    %are the same size), and then computes static and dynamic MH loops
    %   Properties (State):
    %   Methods:
    
    %% Propterties
    properties (SetAccess = protected) %Set access is from class or subclasses
        particles
        R               %interparticle radius vector
        r               %interparticle radius amplitude norm(R)
        U_part          %cell array of particle energies (anon. functions)
        U_total         %total energy for model (anon. function)
        ode_options     %options for ODE solver
        ode_solver      %indicates which ode solver to use
        opt_options     %options for energy min finder        
        flag_nearest    %boolean identifier of whether or not to dipole couple
        ode_solver_list %list of available ode sovlers
        
    end
    properties (Access = public)     %Access is from class or subclasses
        %these are primarily used as temporary variables for the solvers (mainly
        %the LLG solver)
        Plots           %access the visualization class MSPlot_
        count
        progress
        time
        ang
        m
        Heff
        
        r_cutoff        %nearest neighbor cuttoff distance
        static_solver_type     %coupled or segregated 
        seg_tol         %angular tolerance for segmented solver
    end
    
    %% Methods
    methods
        %% Constructor
        function obj = MSModel_(particles)
            %make sure particles is a cell array full of MSparticles_
            switch iscell(particles);
                case 0
                    error('Particles must be stored in a cell array')
                case 1
                    for i = 1:numel(particles)
                        switch isa(particles{i},'MSParticle_');
                            case 0
                                error('Particle %d is not a %s',i,class(MSParticle_));
                            case 1
                                particles{i}.Update;    %update particle
                        end 
                    end
            end
            
            %initialize variables
            obj.particles = particles;
            obj.Particle_data_check;
            obj.count = 1;
            obj.progress = 0;
            obj.r_cutoff = inf;     %default, include all particles
            
            %set ode_properties
            obj.set_ode_options('default');
            %set opt_options
            obj.set_static_solver_options('default');
            %set static solver type / properties
            obj.set_static_solver_type('default');
            obj.seg_tol = 1e-6; %default segmented solver tolerance
            %set dynamic solver type and list available solvers
            obj.ode_solver = 'ode45'; 
            obj.ode_solver_list = {'ode45','ode23','ode113','ode15s'};
            
            %access to built in plots
            obj.Plots = MSPlot_(obj);
            
        end
        
        %% Set ODE options
        function obj = set_ode_options(varargin)
            obj = varargin{1};
            switch nargin > 2 || ~strcmp(varargin{2},'default')
                case 0 %load default optoins
                    obj.ode_options = odeset(...
                        'RelTol',1e-6,...
                        'AbsTol',1e-9,...
                        'NormControl','on',...
                        'InitialStep',1e-14,...
                        'MaxStep',1e-10);
                case 1 %assign input options
                    %separate input into property / value pairs
                    [prop, val] = PropertyValue(varargin(2:end));
                    
                    %get 'solver' if one was entered
                    solv_ind = strcmp(prop,'solver');
                    switch any(solv_ind)
                        case 1
                            %make sure solver is in allowed list
                            switch any(ismember(obj.ode_solver_list, val{solv_ind}))
                                case 1
                                    obj.ode_solver = val{solv_ind};
                                case 0
                                    temp = obj.ode_solver_list;
                                    temp(2,:) = {', '};
                                    temp(end,end) = {' '};
                                    error('%s is not an allowed solver. The currently implemented solvers are:\n    %s',val{solv_ind},[temp{:}])
                            end
                            %remove solver from prop/val pairs
                            prop = prop(~solv_ind);
                            val = val(~solv_ind);
                        case 0
                    end

                    
                    %initialize
                    if isempty(obj.ode_options)
                        obj.ode_options = odeset;
                    end
                    %set options
                    for i = 1:numel(prop)
                        obj.ode_options.(prop{i}) = val{i};
                    end
            end
            
        end
        
        %% Set fmin (optimization) options
        function obj = set_static_solver_options(varargin)
            obj = varargin{1};
            switch nargin > 2 || ~strcmp(varargin{2},'default')
                case 0 %load default optoins
                    obj.opt_options = optimset;
                    obj.opt_options = optimset(...
                        'Display','off',...
                        'TolX',1e-4,...
                        'TolFun',1e-4);
                case 1 %assign input options
                    %separate input into property / value pairs
                    [prop, val] = PropertyValue(varargin(2:end));
                    
                    %initialize
                    if isempty(obj.ode_options)
                        obj.opt_options = optimset;
                    end
                    %set options
                    for i = 1:numel(prop)
                        obj.opt_options.(prop{i}) = val{i};
                    end
            end
            
        end
        
        %% Set static solver type
        function obj = set_static_solver_type(varargin)
            obj = varargin{1};
            switch nargin == 1 || strcmp(varargin{2},'default')
                case 1 %load default optoins
                    obj.static_solver_type = 'segregated';
                case 0 %assign input options
                    switch varargin{2}
                        case 'segregated'
                            obj.static_solver_type = 'segregated';
                        case 'coupled'
                            obj.static_solver_type = 'coupled' ;
                        otherwise
                            error('Available solver types are "segregated" or "coupled".\n"%s" is now an allowed solver type',varargin{2})
                    end
            end
            
        end
        
        %% Calculations (MHstatic, MHdynamic, Heff, LLG)
        function obj = RUN_static(obj)
            obj = interparticle_dist(obj); %update interparticle distance            
            obj = Static_MH(obj);            
        end
        
        function obj = RUN_LLG(obj)
            obj = interparticle_dist(obj); %update interparticle distance
            obj.progress = 0;
            obj.count = 1;
            obj = Dynamic_MH_LLG(obj);
        end
        
    end
    
    methods (Access = private)
        %% checkout particle values
        function obj = Particle_data_check(obj)
            %parse input properties
            particle = obj.particles;
            static_props = particle{1}.SourceFields.static.Property_List;
            dynamic_props = particle{1}.SourceFields.dynamic.Property_List;
            
            %initialize variables
            length_static = numel( particle{1}.SourceFields.static.(static_props{1}));
            length_dynamic = numel( particle{1}.SourceFields.dynamic.(dynamic_props{1}));
            
            %compare data
            for i = 1:numel(particle)
                %static
                for j = 1:numel(static_props)
                    temp = numel( particle{i}.SourceFields.static.(static_props{j}));
                    switch temp == length_static
                        case 0
                            error('The input static fields for particle %d \nneed to have the same size as the others',i)
                    end
                end
                
                %dynamic
                for j = 1:numel(dynamic_props)
                    temp = numel( particle{i}.SourceFields.dynamic.(dynamic_props{j}));
                    switch temp == length_dynamic
                        case 0
                            error('The input dynamic fields for particle %d \nneed to have the same size as the others',i)
                    end
                end
                
            end
        end
        
    end
    
    %% call solvers
    methods        
        
        %% Static Equilibrium Behavior
        function [ model ] = Static_MH( model )
            %STATIC_MH Summary of this function goes here
            %   The function finds the energy associated with the magnetic orientation,
            %   and finds the nearest local minimum / equilibrium point
            
            % Setup
            particles = model.particles;    % get particles
            nparts = numel(particles);      % number of particles
            
            % Assemble Total Field functions
            [H1_tot,H2_tot, H3_tot] = total_field(model,nparts);
            
            % Assembled particle energy functions
            model = particle_energy(model,particles,H1_tot,H2_tot,H3_tot,nparts);
            
            % Assemble total energy function
            model = total_energy(model,particles,H1_tot,H2_tot,H3_tot,nparts);
            
            %solve for energy minima
            switch model.static_solver_type
                case 'segregated'
                    [theta,phi] = segregated_solver(model);
                case 'coupled'
                    [theta,phi] = coupled_solver(model);
                otherwise
                    error('Available solver types are "segregated" or "coupled".\n %s is now an allowed solver type',model.static_solver_type)
            end
            
            % Convert back to direction cosines
            [m1,m2,m3] = sph2cart(phi,pi/2-theta,ones(size(phi)));
            for i = 1:nparts
                particles{i}.State.set_state(m1(:,i),'Stat_m1');
                particles{i}.State.set_state(m2(:,i),'Stat_m2');
                particles{i}.State.set_state(m3(:,i),'Stat_m3');
            end
            
        end%Static_MH_Dipole()
        
        %% LLG / Dynamic Solver
        function [ model ] = Dynamic_MH_LLG( model )
            nparts = numel(model.particles);      %number of particles
            time = model.particles{1}.SourceFields.dynamic.t; %#ok<*PROP> %specify specfic time points
            
            %initial magnetization state (3 x n) and angle state (2 x n)
            m0 = zeros(3,nparts); %(m1,m2,m3)
            ang0 = zeros(2,nparts);%(theta, phi)
            for i = 1:nparts
                m0(1,i) = model.particles{i}.State.m1;                            %m1
                m0(2,i) = model.particles{i}.State.m2;                            %m2
                m0(3,i) = model.particles{i}.State.m3;                            %m3
                ang0(1,i) = atan2(sqrt(m0(1,i)^2 + m0(2,i)^2) , m0(3,i) );  %theta
                ang0(2,i) = atan2(m0(2,i), m0(1,i) );                       %phi
            end
            
            p = ProgressBar(1e3);%progress bar
                       
            % Create anonymous ODE function to be solved
            dang_dt = @(t,ang) LLG_angular_equation(t,ang,model,model.particles, nparts,time,p);
            %solve LLG equation
            
            switch model.ode_solver
                case 'ode45'
                    [t_solved,ang_solved]  = ode45(dang_dt, time, ang0(:), model.ode_options);
                case 'ode23'
                    [t_solved,ang_solved]  = ode23(dang_dt, time, ang0(:), model.ode_options);
                case 'ode113'
                    [t_solved,ang_solved]  = ode113(dang_dt, time, ang0(:), model.ode_options);
                case 'ode15s'
                    [t_solved,ang_solved]  = ode15s(dang_dt, time, ang0(:), model.ode_options);
            end
            
            p.stop; %stop progress bar
            
            %ensure t_solved = time
            switch all(t_solved(:) == time(:));
                case 0
                    warning('ode solver output different times than those input')
            end
            
            %% Update Particles with solution
            for i = 0:nparts-1
                model.particles{i+1}.State.set_state(sin(ang_solved(:,2*i+1)).*cos(ang_solved(:,2*i+2)),'Dyn_m1');
                model.particles{i+1}.State.set_state(sin(ang_solved(:,2*i+1)).*sin(ang_solved(:,2*i+2)),'Dyn_m2');
                model.particles{i+1}.State.set_state(cos(ang_solved(:,2*i+1)),'Dyn_m3');
            end                       
        end %Dynamic_MH_LLG()              
                
    end %Methods
        
end%MSParticle_

%%
%{%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Local Functions: Available to the Class only            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%}


%% Interparticle Distances
function obj = interparticle_dist(obj)
particles = obj.particles;
nparts = numel(particles);
%find distance for particles
R = cell(nparts, nparts);  %X_col - X_row : points from row to col
r = zeros(size(R));          %norm(R)
flag = true(size(R));

for i = 1:nparts
    for j = 1:nparts
        R{i,j} = particles{j}.Properties.Location - particles{i}.Properties.Location;
        r(i,j) = norm(R{i,j});                
        %modify based on nearest neighbor distance
        switch r(i,j) > obj.r_cutoff || r(i,j) == 0        
            case 1
                flag(i,j) = false; %don't dipole include coupling for particles at zero distance
        end        
    end    
end
obj.R = R;
obj.r = r;
obj.flag_nearest = flag;
end

%% Total effective field, including dipole coupling
function [H1_tot,H2_tot, H3_tot] = total_field(model,nparts)
particles = model.particles;
%anonymous functions
MsV = @(i) particles{i}.Properties.Volume .* particles{i}.Properties.Ms;
%find distance for particles / flag indicating if particle is too far awaw
R = model.R;
r = model.r;
flag = model.flag_nearest;

%Dipole field H_demag{i,j} - field from particle i on particle j
H_dipole = cell(size(R));
H1_dipole = cell(size(R));
H2_dipole = cell(size(R));
H3_dipole = cell(size(R));
for i = 1:nparts
    for j = 1:nparts
        switch flag(i,j)
            case 1 %include dipole coupling
                H_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(:).*MsV(i).*...
                    (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    MsV(i).*[M1;M2;M3]./(r(i,j)+eps).^3  );
                H1_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(1).*MsV(i).*...
                    (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    MsV(i).*M1./(r(i,j)+eps).^3  );
                H2_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(2).*MsV(i).*...
                    (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    MsV(i).*M2./(r(i,j)+eps).^3  );
                H3_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(3).*MsV(i).*...
                    (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    MsV(i).*M3./(r(i,j)+eps).^3  );
            case 0 %turn off dipole coupling
                H_dipole{i,j} = @(M1,M2,M3) [0;0;0];
                H1_dipole{i,j} = @(M1,M2,M3) 0;
                H2_dipole{i,j} = @(M1,M2,M3) 0;
                H3_dipole{i,j} = @(M1,M2,M3) 0;
        end
    end
end

% Total field
%calculate total dipole field on particle i, as a function of the
%angles of the other particles

%angles will be an array (1, 2*nparts), where the first nparts numbers
%correspond to theta for each particle, and the second nparts numbers
%correspond to phi

%direction cosines
M1 = @(theta,phi) sin(theta).*cos(phi); %#ok<NASGU>
M2 = @(theta,phi) sin(theta).*sin(phi); %#ok<NASGU>
M3 = @(theta,phi) cos(theta);           %#ok<NASGU>

%demag field  H_dipole{j,i}(M1,M2,M3) Field from j on i
H1_tot = cell(nparts,1);
H2_tot = cell(nparts,1);
H3_tot = cell(nparts,1);

for i = 1:nparts
    %Applied field at time step t
    eval1 = ['H1_tot{',num2str(i),'} = @(t,angles) particles{',num2str(i),'}.SourceFields.static.H1(t)'];
    eval2 = ['H2_tot{',num2str(i),'} = @(t,angles) particles{',num2str(i),'}.SourceFields.static.H2(t)'];
    eval3 = ['H3_tot{',num2str(i),'} = @(t,angles) particles{',num2str(i),'}.SourceFields.static.H3(t)'];
    %add dipole fields
    for j = 1:nparts
        switch flag(i,j)
            case 1
                eval1 = [eval1, ' + ','H1_dipole{',num2str(j),',',num2str(i),'}(M1(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M2(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M3(angles(',num2str(j),'),angles(nparts+',num2str(j),')) )' ];  %#ok<AGROW>
                eval2 = [eval2, ' + ','H2_dipole{',num2str(j),',',num2str(i),'}(M1(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M2(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M3(angles(',num2str(j),'),angles(nparts+',num2str(j),')) )' ];  %#ok<AGROW>
                eval3 = [eval3, ' + ','H2_dipole{',num2str(j),',',num2str(i),'}(M1(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M2(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M3(angles(',num2str(j),'),angles(nparts+',num2str(j),')) )' ];  %#ok<AGROW>
        end        
    end
    %evaluate / create the anonymous function 
    eval([eval1,';'])
    eval([eval2,';'])
    eval([eval3,';'])
end

end

%% Numeric Total effective field, including dipole coupling
function [H_tot] = numeric_total_field(model,m,time,t,nparts)
%similar to "total_field", but it takes the current magnetization state,
%and computes numeric values of all effective fields, including the dipole
%field
%m is a (3 x nparts) matrix  containing the m vector for all particles

particles = model.particles;
%anonymous functions
Vol = @(i) particles{i}.Properties.Volume;
Ms = @(i) particles{i}.Properties.Ms;

%Magnetizations
M1 = @(i) Vol(i)*Ms(i)* m(1,i);
M2 = @(i) Vol(i)*Ms(i)* m(2,i);
M3 = @(i) Vol(i)*Ms(i)* m(3,i);

%find distance for particles / flag indicating if particle is too far awaw
R = model.R;
r = model.r;
flag = model.flag_nearest;

%Dipole field H_dipole(i,j) - field from particle i on particle j
H1_dipole = zeros(size(R));
H2_dipole = H1_dipole;
H3_dipole = H1_dipole;
for i = 1:nparts
    for j = 1:nparts
        switch flag(i,j)
            case 1 %include dipole coupling
                H1_dipole(i,j) = 1/(4*pi) .*(3.*R{i,j}(1).*...
                    (M1(i).*R{i,j}(1) + M2(i).*R{i,j}(2) + M3(i).*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    M1(i)./(r(i,j)+eps).^3  );
                H2_dipole(i,j) = 1/(4*pi) .*(3.*R{i,j}(2).*...
                    (M1(i).*R{i,j}(1) + M2(i).*R{i,j}(2) + M3(i).*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    M2(i)./(r(i,j)+eps).^3  );
                H3_dipole(i,j) = 1/(4*pi) .*(3.*R{i,j}(3).*...
                    (M1(i).*R{i,j}(1) + M2(i).*R{i,j}(2) + M3(i).*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
                    M3(i)./(r(i,j)+eps).^3  );
            case 0 %no dipole coupling
                H1_dipole(i,j) = 0;
                H2_dipole(i,j) = 0;
                H3_dipole(i,j) = 0;
        end
    end
end
%sum all the dipole fields on each particle
% H1_dipole = sum(H1_dipole,1);
% H2_dipole = sum(H2_dipole,1);
% H3_dipole = sum(H3_dipole,1);
H_dipole = [sum(H1_dipole,1);sum(H2_dipole,1);sum(H3_dipole,1)];

% Local Effective Field - interpolated at given time step
H_eff = zeros(3,nparts);
for i = 1:nparts
    f = particles{i}.SourceFields.dynamic;
    fint = @(name,t) interp1(time,f.(name),t);   %when time=t, find the value of name
    H_eff(1,i) = particles{i}.EffectiveFields.H1_total(...
        fint('H1',t),fint('H2',t),fint('H3',t),...
        m(1,i),m(2,i),m(3,i),...
        fint('S1',t),fint('S2',t),fint('S3',t),fint('S4',t),fint('S5',t),fint('S6',t),...
        fint('sigma1',t),fint('sigma2',t),fint('sigma3',t));
    H_eff(2,i) = particles{i}.EffectiveFields.H2_total(...
        fint('H1',t),fint('H2',t),fint('H3',t),...
        m(1,i),m(2,i),m(3,i),...
        fint('S1',t),fint('S2',t),fint('S3',t),fint('S4',t),fint('S5',t),fint('S6',t),...
        fint('sigma1',t),fint('sigma2',t),fint('sigma3',t));
    H_eff(3,i) = particles{i}.EffectiveFields.H3_total(...
        fint('H1',t),fint('H2',t),fint('H3',t),...
        m(1,i),m(2,i),m(3,i),...
        fint('S1',t),fint('S2',t),fint('S3',t),fint('S4',t),fint('S5',t),fint('S6',t),...
        fint('sigma1',t),fint('sigma2',t),fint('sigma3',t));    
end


% Total field
H_tot = H_dipole + H_eff;
% H1_tot = H1_dipole + H_eff(1,:);
% H2_tot = H2_dipole + H_eff(2,:);
% H3_tot = H3_dipole + H_eff(3,:);

end

%% Assemble Particle Energy Functions
function model = particle_energy(model,particles,H1_tot,H2_tot,H3_tot,nparts)
% angles: 1 x 2n array - the first nparts are theta, the second are phi
% time: scalar index indicating which value of the applied fields to use

%direction cosines
M1 = @(theta,phi) sin(theta).*cos(phi);
M2 = @(theta,phi) sin(theta).*sin(phi);
M3 = @(theta,phi) cos(theta);

%initialize function
U_part = cell(nparts,1);
for i = 1:nparts
    %total energy function for each particle
    U_part{i} = @(angles,time) particles{i}.Energy.U_total(...
        H1_tot{i}(time,angles),...
        H2_tot{i}(time,angles),...
        H3_tot{i}(time,angles),...
        M1(angles(i),angles(nparts+i)),...
        M2(angles(i),angles(nparts+i)),...
        M3(angles(i),angles(nparts+i)),...
        particles{i}.SourceFields.static.S1(time),...
        particles{i}.SourceFields.static.S2(time),...
        particles{i}.SourceFields.static.S3(time),...
        particles{i}.SourceFields.static.S4(time),...
        particles{i}.SourceFields.static.S5(time),...
        particles{i}.SourceFields.static.S6(time) );
end
model.U_part = U_part;

end

%% Assemble total energy function (all particles)
function model = total_energy(model, particles,H1_tot,H2_tot,H3_tot,nparts) %#ok<INUSL> (suppress warning saying variables aren't used)

%direction cosines
M1 = @(theta,phi) sin(theta).*cos(phi); %#ok<NASGU>
M2 = @(theta,phi) sin(theta).*sin(phi); %#ok<NASGU>
M3 = @(theta,phi) cos(theta);           %#ok<NASGU>

U_part = model.U_part;

%assemble total energy
eval_str = 'U_part{1}(angles,time)';
for i = 2:nparts
    eval_str = [eval_str, ' + ','U_part{',num2str(i),'}(angles,time)'; ]; %#ok<AGROW>
end

for i = 1:nparts
    switch i==1
        case 1 
            eval_str = ['particles{',num2str(i),'}.Energy.U_total( H1_tot{',num2str(i),'}(time,angles), H2_tot{',num2str(i),'}(time,angles),H3_tot{',num2str(i),'}(time,angles),M1(angles(',num2str(i),'),angles(',num2str(nparts+i),')),M2(angles(',num2str(i),'),angles(',num2str(nparts+i),')),M3(angles(',num2str(i),'),angles(',num2str(nparts+i),')),particles{',num2str(i),'}.SourceFields.static.S1(time),particles{',num2str(i),'}.SourceFields.static.S2(time),particles{',num2str(i),'}.SourceFields.static.S3(time),particles{',num2str(i),'}.SourceFields.static.S4(time),particles{',num2str(i),'}.SourceFields.static.S5(time),particles{',num2str(i),'}.SourceFields.static.S6(time))'];
        case 0 
            eval_str = [eval_str,' + particles{',num2str(i),'}.Energy.U_total( H1_tot{',num2str(i),'}(time,angles), H2_tot{',num2str(i),'}(time,angles),H3_tot{',num2str(i),'}(time,angles),M1(angles(',num2str(i),'),angles(',num2str(nparts+i),')),M2(angles(',num2str(i),'),angles(',num2str(nparts+i),')),M3(angles(',num2str(i),'),angles(',num2str(nparts+i),')),particles{',num2str(i),'}.SourceFields.static.S1(time),particles{',num2str(i),'}.SourceFields.static.S2(time),particles{',num2str(i),'}.SourceFields.static.S3(time),particles{',num2str(i),'}.SourceFields.static.S4(time),particles{',num2str(i),'}.SourceFields.static.S5(time),particles{',num2str(i),'}.SourceFields.static.S6(time))']; %#ok<AGROW>
    end
end
eval_str = ['model.U_total = @(angles,time) ', eval_str,';'];
eval(eval_str); %creates U_tot(angles,time)

end

%% Static Solver - Segregated 
function [theta,phi] = segregated_solver(model)
% The segregated solver assumes there is realtively weak coupling between
% neighboring particles. Using this approach, the local min for each
% particle is found based on the current state of the neighboring
% particles. The local min is then iterated over until the angles converge
% to set values

particles = model.particles;    % Setup
nparts = numel(particles);      % number of particles
npts = numel(particles{1}.SourceFields.static.H1); % number of field points
U_part = model.U_part;          % Particle energy functions @(angles, time)
U_tot = model.U_total;

%initialize output
theta = zeros(npts,nparts);
phi = zeros(npts,nparts);
theta0 = zeros(1,nparts);
phi0 = theta0;

p = ProgressBar(npts);   % setup progress bar                
for i = 1:npts          %loop over all field points, and find local minimum state
    converged = false; 
    count = 1;
    while ~converged
        switch i
            case 1 %get initial orientation using fully coupled solver
                for k = 1:nparts
                    theta0(k) = acos(particles{k}.State.m3);
                    phi0(k) = atan2(particles{k}.State.m2,particles{k}.State.m1);
                end
                angles0 = [theta0, phi0];
                U_tot_i = @(angles) U_tot(angles,i);
                angles = fminsearch(U_tot_i,angles0);
                U_tot_temp0 = U_tot_i(angles);
                converged = true;               
                
            otherwise % segregated solver
                switch count > 1
                    case 0 %get initial points from previous time step
                        theta0 = theta(i-1,:);
                        phi0 = phi(i-1,:);
                    case 1 %update points for current time step
                        theta0 = theta(i,:);
                        phi0 = phi(i,:);
                end
                angles0 = [theta0, phi0];
                temp_theta = zeros(1,nparts);
                temp_phi = temp_theta;
                U_tot_temp = 0;
                for j = 1:nparts    %find local min for each particle
                    angles_fn = @(new_angle) [...
                        angles0(1:j-1), new_angle(1), angles0(j+1:nparts),... %theta
                        angles0(nparts+1:nparts+j-1), new_angle(2), angles0(nparts+j+1:end)... %phi
                        ];                    
                    U_part_j = @(new_angle) U_part{j}(angles_fn(new_angle),i);     %particle energy at current time step
                    temp_angles = fminsearch(U_part_j,...       %Find equilibrium orientation
                        [theta0(j)+(rand(1)-0.5)*1e-6,...       %points from preceeding step are slightly perturbed. early models were getting stuck in highly localized single points for no god damned good reason. This fixed it
                        phi0(j)+(rand(1)-0.5)*1e-6],...
                        model.opt_options);                       
                    temp_theta(j) = temp_angles(1);     %store angles in temporary variable
                    temp_phi(j) = temp_angles(2);                  
                    U_tot_temp = U_tot_temp + U_part_j(temp_angles);
                end
                % check if results have converged
                angles = [temp_theta temp_phi];                
                check_convergence = abs(U_tot_temp - U_tot_temp0); %checking energy, instead of angles, is robust to cases where changing an angle doesn't matter (e.g., phi can be any angle when theta is 0)
                U_tot_temp0 = U_tot_temp;
                switch (count > 1) && (check_convergence < model.seg_tol || count > 100 )
                    case 0
                        converged = false;
                    case 1
                        converged = true;                                                               
                end                
        end
        %store angles
        theta(i,:) = mod(angles(1:nparts),2*pi);
        phi(i,:) = mod(angles(nparts+1:end),2*pi);                
        count = count + 1; 
    end
    p.progress; %update progress bar
end
p.stop;             %stop progress bar

end

%% Static Solver - Coupled
function [theta,phi] = coupled_solver(model)
%All angles are optimized at once. It takes longer for a single minimum to
%be found, as you go from minimizing a function of two variables
%(theta,phi), to minimizing one with 2^n_particles variables. While it is
%slower, it is deinitely more robust.

particles = model.particles;    % Setup
nparts = numel(particles);      % number of particles
npts = numel(particles{1}.SourceFields.static.H1); % number of field points
U_part = model.U_part;          % Particle energy functions @(angles, time)
U_tot = model.U_total;            % Total system energy @(angles, time)

%initialize output
theta = zeros(npts,nparts);
phi = zeros(npts,nparts);
theta0 = zeros(1,nparts);
phi0 = theta0;

p = ProgressBar(npts);  % setup progress bar
% Find equilibrium orientations
for i = 1:npts %loop over all field points, and find local minimum state
    for j = 1:nparts        
        %get initial orientation of each particle
        switch i == 1
            case 1 %first point getst the input state                
                theta0(j) = acos(particles{j}.State.m3);
                phi0(j) = atan2(particles{j}.State.m2,particles{j}.State.m1);
            case 0 %subsequent points get values from previous time step
                theta0(j) = angles(j);
                phi0(j) = angles(j+nparts);
        end
    end
    angles0 = [theta0, phi0];               %store initial angles
    U_tot_i = @(angles) U_tot(angles,i);    %store energy at given time   
    [angles] = fminsearch(U_tot_i,angles0+(rand(1)-0.5)*1e-6,model.opt_options);   %Find equilibrium orientation
    
    theta(i,:) = angles(1:nparts);          %store equilibrium angles
    phi(i,:) = angles(nparts+1:end);
    
    p.progress;         %update progress bar
end
p.stop;                 %stop progress bar

end

%% LLG Equation - using angular coordinates
function  [dangdt,m,H_tot] = LLG_angular_equation(t,ang_in,model,particles,nparts,time,pb)
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
[H_tot] = numeric_total_field(model,m,time,t,nparts);

%% Compute dm/dt
%store factors for each particle
F1 = zeros(3,nparts);
F2 = F1;
for i = 1:nparts
    F1(1:3,i) = particles{i}.Properties.Factor_1 ;
    F2(1:3,i) = particles{i}.Properties.Factor_2 ;
end
%Vector form
dmdt = - F1 .* cross(m,H_tot) - F2 .* cross(m,cross(m,H_tot)) ;

%% Convert to dang/dt
%convert m_dot to theta_dot and phi_dot (use chain rule to get m_dot in
%terms of theta_dot and phi_dot, then invert the transformation matrix to
%get these formulas

dangdt = zeros(2,nparts);
%theta_dot
dangdt(1,:) = sec(ang(1,:)).*cos(ang(2,:)).*dmdt(1,:) + ...
    sec(ang(1,:)).*sin(ang(2,:)).*dmdt(2,:);
%phi_dot
if any(ang(1,:) == 0) %avoid numerical inconsistencies when theta = 0 in csc terms
   ang(1,:) = ang(1,:) + (rand(size(ang(1,:)))-.5)*eps; 
end
dangdt(2,:) = -csc(ang(1,:)).*sin(ang(2,:)).*dmdt(1,:) + ...
    csc(ang(1,:)).*cos(ang(2,:)).*dmdt(2,:);


%% update model
model.time{model.count} = t; 
model.m{model.count} = m;
model.ang{model.count} = ang;
model.Heff{model.count} = H_tot;
model.count = model.count + 1;

%% output variables
dangdt = dangdt(:);
m = m(:);
H_tot = H_tot(:);

%% Update progress bar
progress = floor((t-time(1)) / (time(end) - time(1))*1e3);
dif = progress - model.progress;
switch dif > 0
    case 1
        for i = 1:dif
            pb.progress;
            model.progress = model.progress + 1;
        end
    case 0
end

end
