%% M?todos matem?ticos
%% S?ries de Fourier

clear all;
close all;

%% Exemplo 1 - Fun??o onda quadrada
T=2*pi;
t1=-pi:0.01:pi;
n1=length(t1);

for j=1:n1
    if t1(j)<0
        f1(j)=-1;
    else
        f1(j)=1;
    end
end

t=-3*pi:0.01:3*pi;
f=[f1(1:(n1-1)) f1 f1(2:n1)];


figure(1);
plot(t,f,'LineWidth',2);
grid on;

N=20;
fs=zeros(1,length(t));

%% FFT
%FFT eh de 0 a T, nesse caso deve ser de 0 a 2*pi.
t2=0:0.01:(2*pi);
n2=length(t2);
for j=1:n2
    if t2(j)<pi
        f2(j)=1;
    else
        f2(j)=-1;
    end
end

cp_dft=fft(f2)/n2;      
cp=cp_dft(2:N+1);

for p=1:N;
    
%     if mod(p,2) %Odd
%         cp(p)=-i*2/(p*pi);
%     else %Even
%         cp(p) = 0;
%     end
    
    fs_p=cp(p)*exp(i*p*2*pi*t/T) + cp(p)'*exp(-i*p*2*pi*t/T)
    fs=fs+fs_p;
    %figure(2);
    %plot(t,fs_p);
    %hold on;
    %grid on;
end

figure(2);

subplot(1,2,1);
plot(t,f,'LineWidth',2);
hold on
plot(t,fs,'r');
hold off;
grid on;

subplot(1,2,2);
stem([0:1:N],[0 abs(cp)],'o','filled')
xlabel('Frequency (Hz)')
ylabel('Modulus (\mid\itc_n\rm\mid)')

figure(3);

subplot(1,2,1);
stem([0:1:N],[0 abs(cp)],'o','filled')
xlabel('Frequencia (rad/s)')
ylabel('Modulo - \mid\itc_p\rm\mid ')

subplot(1,2,2);
stem([0:1:N],[0 angle(cp)],'o','filled')
xlabel('Frequencia (rad/s)')
ylabel('Fase - arg(\itc_p\rm)')



%% Sistema
m=10;
k=250;
xi=0.01;

wn=sqrt(k/m);
w=0:0.1:20;

H=1/m./(wn^2-w.^2+i*2*xi*w*wn);

figure(4);

subplot(1,2,1);
semilogy(w,abs(H),'LineWidth',2);
xlabel('Frequency (Hz)')
ylabel('Modulo de H')
grid on;

subplot(1,2,2);
plot(w,angle(H),'LineWidth',2);
xlabel('Frequency (Hz)')
ylabel('Modulo de H')
grid on;


%Resposta
xs=zeros(1,length(t));
xs2=zeros(1,length(t));
for p=1:N;
    
%     if mod(p,2) %Odd
%         cp(p)=-i*2/(p*pi);
%     else %Even
%         cp(p) = 0;
%     end
    Hp(p)=1/m./(wn^2-p.^2+i*2*xi*p*wn);
    xs_p=cp(p)*Hp(p)*exp(i*p*2*pi*t/T) + cp(p)'*Hp(p)'*exp(-i*p*2*pi*t/T)
    xs=xs+xs_p;
    xs_p2=2*abs(cp(p))*abs(Hp(p))*cos(p*2*pi*t/T+angle(cp(p))+angle(Hp(p)));    
    xs2=xs2+xs_p2;
    %figure(2);
    %plot(t,fs_p);
    %hold on;
    %grid on;
end

figure(2);

subplot(1,2,1);
plot(t,f,'LineWidth',2);
hold on
plot(t,fs,'r');
hold off;
grid on;

subplot(1,2,2);
plot(t,xs,'LineWidth',2);
hold on
plot(t,xs2,'r');
hold off;
grid on;



figure(3);

subplot(1,2,1);
stem([0:1:N],[0 abs(cp)],'o','filled')
xlabel('Frequency (Hz)')
ylabel('Modulus (\mid\itc_n\rm\mid)')

subplot(1,2,2);
stem([0:1:N],[0 abs(cp.*Hp)],'o','filled')
xlabel('Frequency (Hz)')
ylabel('Modulus (\mid\itc_n\rm\mid)')


