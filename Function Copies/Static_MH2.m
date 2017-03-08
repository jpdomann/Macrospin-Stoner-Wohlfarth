function [ particles ] = Static_MH2( particles )
%STATIC_MH Summary of this function goes here
%   The function finds the energy associated with the magnetic orientation,
%   and finds the nearest local minimum / equilibrium point

%% Setup
%number of particles
nparts = numel(particles);

%number of field points
npts = numel(particles{1}.SourceFields.static.H1);

%initialize output
theta = zeros(npts,nparts);
phi = zeros(npts,nparts);


%% Assemble Total Field functions
[H1_tot,H2_tot, H3_tot] = total_field(particles,nparts);

%% Assemble total energy function
U_tot = total_energy(particles,H1_tot,H2_tot,H3_tot,nparts);

%% Find equilibrium orientations
% setup progress bar
p = ProgressBar(npts);

%loop over all field points, and find local minimum state
theta0 = zeros(1,nparts);
phi0 = theta0;
for i = 1:npts
    for j = 1:nparts
        switch i == 1
            case 1
                %initial orientation for each particle
                theta0(j) = acos(particles{j}.State.m3);
                phi0(j) = atan2(particles{j}.State.m2,particles{j}.State.m1);
            case 0
                theta0(j) = angles(j);
                phi0(j) = angles(j+nparts);
        end
    end
    angles0 = [theta0, phi0];               %store initial angles
    U_tot_i = @(angles) U_tot(angles,i);    %store energy at given time
    angles = fminsearch(U_tot_i,angles0);   %Find equilibrium orientation
    theta(i,:) = angles(1:nparts);          %store equilibrium angles
    phi(i,:) = angles(nparts+1:end);
    
    %update progress bar
    p.progress;
    
end
%stop progress bar
p.stop;

%% Convert back to direction cosines
[m1,m2,m3] = sph2cart(phi,pi/2-theta,ones(size(phi)));

for i = 1:nparts
    particles{i}.State.set_state(m1(:,i),'Stat_m1');
    particles{i}.State.set_state(m2(:,i),'Stat_m2');
    particles{i}.State.set_state(m3(:,i),'Stat_m3');
end

%% Plot Data
% %create a spherical grid of possble m orientations
% nmesh = 100;
% [THETA, PHI] = meshgrid(linspace(0,pi,nmesh/2),linspace(-pi,pi,100));
% % R = ones(size(THETA));
% % [M1,M2,M3] = sph2cart(PHI,pi/2-THETA,R);
%
% %setup debugging plot
% figure(1)
% clf
% subplot(1,2,1)
% h1 = surf(THETA*180/pi,PHI*180/pi,ones(size(THETA))); hold all
% h1.EdgeAlpha = 0.1; h1.FaceColor = 'interp';
% h2 = plot3(0,0,0,'r*');
% xlabel('\Theta')
% ylabel('\Phi')
% zlabel('Energy (J/m^3')
% xlim([0 180])
% ylim([-180 180])
% view(2)
% grid on
%
% subplot(1,2,2)
% h3 = plot3(0,0,0,'bo'); hold all
% h4 = plot3(0,0,0,'r:');
% xlabel('m1')
% ylabel('m2')
% zlabel('m3')
% xlim([-1 1])
% ylim([-1 1])
% zlim([-1 1])
% grid on
%
% for j = 1:nparts
%     for i = 1:npts
%         Uplot = @(theta,phi) particles{j}.Energy.U_total(...
%             particles{j}.SourceFields.static.H1(i),...
%             particles{j}.SourceFields.static.H2(i),...
%             particles{j}.SourceFields.static.H3(i),...
%             M1(theta,phi),M2(theta,phi),M3(theta,phi),...
%             particles{j}.SourceFields.static.S1(i),...
%             particles{j}.SourceFields.static.S2(i),...
%             particles{j}.SourceFields.static.S3(i),...
%             particles{j}.SourceFields.static.S4(i),...
%             particles{j}.SourceFields.static.S5(i),...
%             particles{j}.SourceFields.static.S6(i) );
%         h1.ZData = (Uplot(THETA,PHI)-min(min(Uplot(THETA,PHI))))/...
%             max(max(abs(Uplot(THETA,PHI)-min(min(Uplot(THETA,PHI))) )));
%         h2.XData = theta(i,j)*180/pi;
%         h2.YData = phi(i,j)*180/pi;
%         h2.ZData = (Uplot(theta(i,j),phi(i,j))-min(min(Uplot(theta(i,j),phi(i,j)))))/...
%             (max(max(abs(Uplot(theta(i,j),phi(i,j))-min(min(Uplot(theta(i,j),phi(i,j)))) )))+eps);
%         h3.XData = m1(i,j);
%         h3.YData = m2(i,j);
%         h3.ZData = m3(i,j);
%         h4.XData = m1(1:i,j);
%         h4.YData = m2(1:i,j);
%         h4.ZData = m3(1:i,j);
%
%         drawnow
%     end
% end

end

function [H1_tot,H2_tot, H3_tot] = total_field(particles,nparts)
%anonymous functions
Vol = @(i) particles{i}.Properties.Volume;
%find distance for particles
R = cell(nparts, nparts-1);  %X_col - X_row : points from row to col
r = zeros(size(R));          %norm(R)

%Demag field H_demag{i,j} - field from particle i on particle j
H_dipole = cell(size(R));
H1_dipole = cell(size(R));
H2_dipole = cell(size(R));
H3_dipole = cell(size(R));
for i = 1:nparts
    for j = [1:i-1,i+1:nparts] % upper / lower triangular regions
        R{i,j} = particles{j}.Properties.Location - particles{i}.Properties.Location;
        r(i,j) = norm(R{i,j});
        H_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(:).*Vol(i).*...
            (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
            Vol(i).*[M1;M2;M3]./(r(i,j)+eps).^3  );
        H1_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(1).*Vol(i).*...
            (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
            Vol(i).*M1./(r(i,j)+eps).^3  );
        H2_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(2).*Vol(i).*...
            (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
            Vol(i).*M2./(r(i,j)+eps).^3  );
        H3_dipole{i,j} = @(M1,M2,M3) 1/(4*pi) .*(3.*R{i,j}(3).*Vol(i).*...
            (M1.*R{i,j}(1) + M2.*R{i,j}(2) + M3.*R{i,j}(3) ) ./ (r(i,j)+eps).^5 - ...
            Vol(i).*M3./(r(i,j)+eps).^3  );
    end
    for j = i %principal diagonal
        H_dipole{i,j} = @(M1,M2,M3) [0;0;0];
        H1_dipole{i,j} = @(M1,M2,M3) 0;
        H2_dipole{i,j} = @(M1,M2,M3) 0;
        H3_dipole{i,j} = @(M1,M2,M3) 0;
    end
end


%% Total field
%calculate total dipole field on particle i, as a function of the
%angles of the other particles

%angles will be an array (1, 2*nparts), where the first nparts numbers
%correspond to theta for each particle, and the second nparts numbers
%correspond to phi

%direction cosines
M1 = @(theta,phi) sin(theta).*cos(phi);
M2 = @(theta,phi) sin(theta).*sin(phi);
M3 = @(theta,phi) cos(theta);

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
        eval1 = [eval1, ' + ','H1_dipole{',num2str(j),',',num2str(i),'}(M1(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M2(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M3(angles(',num2str(j),'),angles(nparts+',num2str(j),')) )' ];
        eval2 = [eval2, ' + ','H2_dipole{',num2str(j),',',num2str(i),'}(M1(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M2(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M3(angles(',num2str(j),'),angles(nparts+',num2str(j),')) )' ];
        eval3 = [eval3, ' + ','H2_dipole{',num2str(j),',',num2str(i),'}(M1(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M2(angles(',num2str(j),'),angles(nparts+',num2str(j),')), M3(angles(',num2str(j),'),angles(nparts+',num2str(j),')) )' ];
    end
    %evaluate
    eval([eval1,';'])
    eval([eval2,';'])
    eval([eval3,';'])
end
end

function U_tot = total_energy(particles,H1_tot,H2_tot,H3_tot,nparts)
%% Assemble total energy function

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
%assemble total energy
eval_str = 'U_tot = @(angles,time) U_part{1}(angles,time)';
for i = 2:nparts
    eval_str = [eval_str, ' + ','U_part{',num2str(i),'}(angles,time)'; ];
end
eval(eval_str); %creates U_tot(angles,time)

end