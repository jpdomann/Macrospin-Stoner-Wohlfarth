function [ particles ] = Static_MH( particles )
%STATIC_MH Summary of this function goes here
%   The function finds the energy associated with the magnetic orientation,
%   and finds the nearest local minimum / equilibrium point

%number of particles
nparts = numel(particles);

%number of field points
npts = numel(particles{1}.SourceFields.static.H1);

%direction cosines
M1 = @(theta,phi) sin(theta).*cos(phi);
M2 = @(theta,phi) sin(theta).*sin(phi);
M3 = @(theta,phi) cos(theta);

%initialize output
theta = zeros(npts,nparts);
phi = zeros(npts,nparts);

% setup progress bar
p = ProgressBar(npts*nparts);

%loop over all field points, and find local minimum state

%Update so all particles are computed at once for a given field point (this
%is needed when particles become coupled together)
for i = 1:npts
    
    for j = 1:nparts
        %total energy function for each particle
        Usym = @(angles) particles{j}.Energy.U_total(...
            particles{j}.SourceFields.static.H1(i),...
            particles{j}.SourceFields.static.H2(i),...
            particles{j}.SourceFields.static.H3(i),...
            M1(angles(1),angles(2)),...
            M2(angles(1),angles(2)),...
            M3(angles(1),angles(2)),...
            particles{j}.SourceFields.static.S1(i),...
            particles{j}.SourceFields.static.S2(i),...
            particles{j}.SourceFields.static.S3(i),...
            particles{j}.SourceFields.static.S4(i),...
            particles{j}.SourceFields.static.S5(i),...
            particles{j}.SourceFields.static.S6(i) );
        
        %initial orientation for each particle
        theta0 = acos(particles{j}.State.m3);
        phi0 = atan2(particles{j}.State.m2,particles{j}.State.m1);
        
        %Find equilibrium orientation
        angles = fminsearch(Usym,[theta0,phi0]);
        theta(i,j) = angles(1);
        phi(i,j) = angles(2);
        
        %update progress bar
        p.progress;
    end
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

% %% Plot Data
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

