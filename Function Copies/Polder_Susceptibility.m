%% LLG Polder Susceptibility
clear
clc

%% Constants
mu0 = 4*pi*1e-7;    %[H/m] vacuum permeability
gamma = 1.760859e11;  %[1/(s*T)] Electron Gyromagnetic ratio              
alpha = 1e-3;
Ms = 1.8 / mu0 ;

%% functions
w0 = @(H) gamma*mu0*H;
wM = w0(Ms);
wr = @(H) w0(H);
% chi1 = @(w,H) wM .* w0(H).*( w0(H).^2 - w.^2.*(1-alpha^2) ) ./...
%     ( (w0(H).^2 - w.^2.*(1+alpha^2) ).^2 + 4.* w.^2 .* w0(H).^2 .* alpha.^2  ) ;
chi1 = @(w,H) wM .* w0(H).*( w0(H).^2 - w.^2.*(1-alpha^2) ) ./...
    ( (w0(H).^2 - w.^2.*(1+alpha^2) ).^2 + 4.* alpha.^2 .* w.^2 .* w0(H).^2  ) ;
chi2 = @(w,H) alpha.* wM .* w .*( w0(H).^2 + w.^2.*(1+alpha^2) ) ./...
    ( (w0(H).^2 - w.^2.*(1+alpha^2) ).^2 + 4.* alpha.^2 .* w.^2 .* w0(H).^2  ) ;

fL = @(H) gamma / (2*pi) * sqrt(mu0.*(H + Ms) .*(mu0.*(H + Ms) ) );

%% Mesh
[f,H] = meshgrid(logspace(8.5,9.5,1e3), linspace(1e4,Ms/20,1e3));

%% Plot
fig = figure(1);
clf
ax = gca;

h = mesh(f,H,chi1(f*2*pi,H));
hold on
% h2 = plot3(fL(unique(H)),unique(H), chi1(2*pi*fL(unique(H)),unique(H) ),'r*');
h2 = plot3(w0(unique(H))/(2*pi),unique(H), chi1(w0(unique(H)),unique(H) ),'r*');
ax.XScale = 'log';
xlabel('f (Hz)')
ylabel('H (A/m)')
zlabel('\chi')

% Line plot
f2 = logspace(8.44,8.55,1e4);
H2 = 9e3;

figure(2)
clf
ax = gca;
h4 = plot(f2,chi2(2*pi*f2,H2),'b');
hold on
h3 = plot(f2,chi1(2*pi*f2,H2),'r');
% ax.XScale = 'log';
ax.FontSize = 12;
ax.FontWeight = 'bold';
h3.LineWidth = 2;
h4.LineWidth = 2;
xlabel('f (Hz)')
ylabel('\mu_r')
legend('Im[\mu_r]','Re[\mu_r]')
grid on
