%% M?todos matem?ticos
%% S?ries de Fourier

clear all;

%% Exemplo 1 - Fun??o onda quadrada

t=-pi:0.01:pi;

for j=1:length(t)
    if t(j)<0
        f(j)=-1;
    else
        f(j)=1;
    end
end

N=100;
fs=zeros(1,length(t));

for n=1:N;
    fs_p=4/(pi*(2*n-1))*sin((2*n-1)*t);
    fs=fs+fs_p;
    figure(1);
    plot(t,fs_p);
    hold on;
    grid on;
end

figure(2);
plot(t,f,'LineWidth',2);
hold on
plot(t,fs,'r');
hold off;
grid on;

% Incluindo uma diferenca de fase
% Fase qualquer
% fs=0;
% for n=1:N;
%     fs_p=4/(pi*(2*n-1))*sin((2*n-1)*t + pi/2);
%     fs=fs+fs_p;
%     figure(3);
%     plot(t,fs_p);
%     hold on;
%     grid on;
% end
% 
% %Deslocamento lateral
% fs=0;
% for n=1:N;
%     fs_p=4/(pi*(2*n-1))*sin((2*n-1)*t + (2*n-1)*pi/2);
%     fs=fs+fs_p;
% end
% 
% 
% figure (4);
% plot(t,f,'LineWidth',2);
% hold on
% plot(t,fs,'r');
% hold off;
% grid on;
% 


% % Exemplo 2 - Fun??o f(x)=x
% 
% x=0:0.001:1;
% 
% f=x;
% 
% N=32;
% fs=ones(1,length(x))/2;
% 
% figure(1);
% plot(x,fs);
% hold on;
% grid on; 
% 
% for j=1:N;
%     fs_p=-1/pi*sin(2*j*pi*x)/j;
%     fs=fs+fs_p;
%     figure(1);
%     plot(x,fs_p);
%     hold on;
%     grid on;    
% end
% 
% figure(2);
% plot(x,f,'LineWidth',2);
% hold on
% plot(x,fs,'r');
% hold off;
% grid on;


