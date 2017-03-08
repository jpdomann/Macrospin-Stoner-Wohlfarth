clear
clc

%% Compute
sp = 2; 

shapes = {'Rectangle','Ellipse'};
shape = shapes{sp};

a=1;
% beta = linspace(0.1, 1, 10);     %b/a points
% gamma = linspace(0.1, 1, 10);    %c/a points
beta = logspace(-4, 0, 50);     %b/a points
gamma = logspace(-4, 0, 50);    %c/a points

[Beta,Gamma] = meshgrid(beta,gamma);
b = a*Beta;
c = a*Gamma;

n = numel(Beta);
%prolate
prolate = zeros(numel(beta),3);
prolate_ba = zeros(size(beta));
prolate_ca = prolate_ba;
for i = 1:size(Beta,2)
    prolate(i,:)=DemagFactor(shape,[a b(1,i) b(1,i)]);
    prolate_ba(i) = b(1,i)/a;
    prolate_ca(i) = b(1,i)/a;
end

%oblate
oblate = zeros(numel(gamma),3);
oblate_ba = zeros(size(gamma));
oblate_ca = prolate_ba;
for i = 1:size(Beta,1)
    oblate(i,:)=DemagFactor(shape,[a a c(i,1)]);
    oblate_ba(i) = a/a;
    oblate_ca(i) = c(i)/a;
end

%general
d = zeros([size(Beta),3]);
for i = 1:size(Beta,1)
    for j = 1:size(Beta,2)
        temp=DemagFactor(shape,[a b(i,j) c(i,j)]);
        d(i,j,1)=temp(1);
        d(i,j,2)=temp(2);
        d(i,j,3)=temp(3);
    end
end

%Sort
N1 = zeros(size(Beta));
N2 = N1; N3 = N1;
for i = 1:size(Beta,1)
    for j = 1:size(Beta,2)
        N = d(i,j,:); N=sort(N(:),'ascend');
        N1(i,j) = N(1);
        N2(i,j) = N(2);
        N3(i,j) = N(3);
    end
end

%% plot

for i = 1:3
    figure(i)    
    subplot(1,2,sp)
    cla
    
    plot(oblate_ca,oblate(:,i),prolate_ca,prolate(:,i),'LineWidth',2)
    hold all
    for j = 1:numel(beta)
        plot(Gamma(Beta == beta(j)),eval(['N',num2str(i),'(Beta == beta(',num2str(j),'))']))
    end
    grid on
end
title(shape)

%% Plot combined
color = {'r','b','g'};
figure(4)
subplot(1,2,sp)
cla

for i = 1:3    
    plot(oblate_ca,oblate(:,i),prolate_ca,prolate(:,i),'LineWidth',2)
    hold all
    for j = 1:numel(beta)
        plot(Gamma(Beta == beta(j)),eval(['N',num2str(i),'(Beta == beta(',num2str(j),'))']),'Color',color{i}  )
    end
    grid on
end
title(shape)

%% Surface 
color = {'r','b','g'};
figure(5)
subplot(1,2,sp)
cla

fa = .4;
mesh(Beta,Gamma,N1,'FaceAlpha',fa,'EdgeAlpha',fa,'FaceColor','b')
hold all
mesh(Beta,Gamma,N2,'FaceAlpha',fa,'EdgeAlpha',fa,'FaceColor','g')
mesh(Beta,Gamma,N3,'FaceAlpha',fa,'EdgeAlpha',fa,'FaceColor','r')
for j = 1:3
    plot3(oblate_ba,oblate_ca,oblate(:,j),'LineWidth',2,'Color',color{j})
    plot3(prolate_ba,prolate_ca,prolate(:,j),'LineWidth',2,'Color',color{j})
end
grid on
legend('N1','N2','N3')
title(shape)
xlabel('B/A')
ylabel('C/A')
zlabel('N')

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
