%%% Sistemas Discretos - Sistema 3 GL
%
% v1 - 17/05/2009
% v2 - 27/06/2013
%
% Resumo:
% - Calcula freq. naturais e formas modais
% - Anima os modos
% - Anima a resposta total
% - Calcula a resposta de uma das massas

clear all;

N=3;

%% Massas
m1=100;      
m2=m1;
m3=m1;

%% Rigidezas
k1=10000;      
k2=k1;
k3=k1;

%% Montando matrizes de massa e rigidez
K=[k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3];
Mv=[m1 m2 m3];
M=diag(Mv);

%% Extra??o dos autovalores e autovetores
[V,W]=eig(K,M);


%% Freq. Naturais e
fi=diag(W.^0.5)/(2*pi);

% figure(1);
% plot(1:length(fi),fi,'s-');
% xlabel('Frequencia natural'); 
% ylabel('Frequencia [Hz]');
% % axis([0 3 -1.2 1.2]);
% grid on



%% Formas modais
figure(2);
set(gca,'OuterPosition',[6.3543e-01  -1.3311e-01   2.7454e-01   1.1932e+00]);
subplot(1,3,1);
plot(0:N,[0; V(:,1)],'-o');
title('Primeiro Modo','Interpreter','latex','FontSize',16);
xlabel('Grau de Liberdade','Interpreter','latex','FontSize',16); 
ylabel('Amplitude','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14);
subplot(1,3,2);
plot(0:N,[0; V(:,2)],'-o');
title('Segundo Modo','Interpreter','latex','FontSize',16);
xlabel('Grau de Liberdade','Interpreter','latex','FontSize',16); 
ylabel('Amplitude','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14);
subplot(1,3,3);
plot(0:N,[0; V(:,3)],'-o');
title('Terceiro Modo','Interpreter','latex','FontSize',16);
xlabel('Grau de Liberdade','Interpreter','latex','FontSize',16); 
ylabel('Amplitude','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14);
print -deps fig_modos


%% Anima??o das Formas modais

mod=1;
figure(3);
for k = 1:3200
    modo=[0 V(1,mod) V(2,mod) V(3,mod)]*sin(1*pi*(k-1)/16);
    modo1=[0 1 2 3]+modo*4;
	plot(modo1,zeros(1,N+1),'-s','markersize',30,'markerfacecolor','b');
	axis([0 4 -1 1])
	mov(k) = getframe;
end
movie(mov,30)


%% Anima??o da resposta do sistema
% 
% t=0:0.01:20;
% peso=[1 1 1]
% figure(4);
% for k = 1:length(t)
%     resp=[0 1 2 3];
%     for m=1:N
%         modo=peso(m)*[0 V(1,m) V(2,m) V(3,m)]*sin(2*pi*fi(m)*t(k));
%         resp=resp+modo*2;
%     end    
% 	plot(resp,zeros(1,N+1),'-s','MarkerSize',30,'MarkerFaceColor','b');
% 	axis([0 4 -1 1])
% 	Mov(k) = getframe;
% end
% 
% movie(Mov,3)

%% Soma resposta modos - Massa 1
% 
% t=0:0.001:10;
% for k=1:length(t)
%     resp_modo1(k)=V(1,1)*sin(2*pi*fi(1)*t(k));
%     resp_modo2(k)=V(1,2)*sin(2*pi*fi(2)*t(k));
%     resp_modo3(k)=V(1,3)*sin(2*pi*fi(3)*t(k));
%     resp_tot(k)=resp_modo1(k)+resp_modo2(k)+resp_modo3(k);
% end
% figure(5);
% plot(t,resp_modo1,'-b');
% hold on
% plot(t,resp_modo2,'-m');
% plot(t,resp_modo3,'-k');
% plot(t,resp_tot,'-r','LineWidth',2);
% xlabel('Tempo [s]'); 
% ylabel('Deslocamento [m]');
% %axis([0 3 -1.2 1.2]);
% grid on
% hold off
% 


