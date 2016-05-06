%% M?todos matem?ticos
%% S?ries de Fourier

clear all;
close all;

%% Excitacao qualquer
t0=1;
F0=1;
t=0:0.01:5;
N=length(t);
for j=1:N
    if t(j)<t0
        f(j)=F0*(1-t(j)/t0);
    else
        f(j)=0;
    end
end


%% Sistema
m=10;
k=5000;
xi=0.03;

wn=sqrt(k/m);
wd=wn*sqrt(1-xi^2);

h=1/(m*wd)*exp(-xi*wn*t).*sin(wd*t);

figure(1);

subplot(1,2,1);
plot(t,f,'LineWidth',2);
xlabel('Tempo [s]');
ylabel('Forca');
grid on;

subplot(1,2,2);
plot(t,h,'LineWidth',2);
xlabel('Tempo [s]')
ylabel('Forca')
grid on;


%Resposta
t2=0:0.01:20;
N2=length(t2);
f2=[f zeros(1,N2-N)];
h2=[h zeros(1,N2-N)];

x=conv(f2,h2);

% figure(2);
% plot(t2(1:N),x(1:N),'LineWidth',2);
% 
% grid on;

pic=figure(3);
plot(t2(1:N),x(1:N),'k','LineWidth',2);
grid on;
xlabel(['$t$ [s]'],'interpreter','latex','FontSize',14)
ylabel(['$x(t)$'],'interpreter','latex','FontSize',14)
axis([0 5 -0.04 0.04]);
hold off;
set(pic,'PaperPositionMode','auto')
%print -deps test_fig



