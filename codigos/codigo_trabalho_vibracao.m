% CODIGO TRABALHO DE VIBRACOES
clc; clear('all'); close all;

numero_graus_liberdade = 4;
% Declaracao das propriedades dos elementos:
% massas:
m_1 = 40; % unidade em Kg
m_2 = m_1;
m_3 = 25; % unidade em Kg
m_4 = m_3;
% molas ou rigidez:
k_1 = 8e6; % unidade em N/m
k_3 = k_1;
k_2 = 2e7; % unidade em N/m
k_5 = k_2;
k_4 = 6e6; % unidade em N/m
k_6 = k_4;
k_7 = 4e6; % unidade em N/m
k_8 = k_7;
% amortecimentos:
c_1 = 360; % unidade em N.s/m
c_3 = c_1;
c_4 = 120; % unidade em N.s/m
c_6 = c_4;
% constantes de amortecimento proporcional:
alpha_proporcional = 0.00001; beta_proporcional = 0.00002;
%alpha_proporcional = 0; beta_proporcional = 0;

% Construindo a matriz de massas
M = [m_1 0 0 0; 0 m_2 0 0; 0 0 m_3 0; 0 0 0 m_4];
% Construindo a matriz de rigidez
K(1,1:4) = [(k_1 + k_2 + k_4), -k_2, -k_4, 0];
K(2,1:4) = [-k_2, (k_3 + k_2 + k_6), 0, -k_6];
K(3,1:4) = [-k_4, 0, (k_5 + k_4 + k_7), -k_5];
K(4,1:4) = [0, -k_6, -k_5, (k_8 + k_6 + k_5)];

% Aplicando amortecimento proporcional para encontrar os outros amortecimentos
c_2 = alpha_proporcional*M(1,1) + beta_proporcional*K(1,1) - c_1 - c_4;
c_5 = alpha_proporcional*M(3,4) + beta_proporcional*K(3,4);
c_7 = alpha_proporcional*M(3,3) + beta_proporcional*K(3,3) - c_5 - c_4;
c_8 = alpha_proporcional*M(4,4) + beta_proporcional*K(4,4) - c_6 - c_5;

% Construindo a matriz de amortecimentos
C(1,1:4) = [(c_1 + c_2 + c_4), -c_2, -c_4, 0];
C(2,1:4) = [-c_2, (c_3 + c_2 + c_6), 0, -c_6];
C(3,1:4) = [-c_4, 0, (c_5 + c_4 + c_7), -c_5];
C(4,1:4) = [0, -c_6, -c_5, (c_8 + c_6 + c_5)];

% Extracao dos autovalores (frequencias naturais) e autovetores (modos de vibracao)
[modos_vibracao, omega_n_quadrado] = eig(K,M);
omega_n = omega_n_quadrado.^0.5;
frequencias_naturais = diag(omega_n)/(2*pi);

% Calculo dos fatores de amortecimento modais
fatores_amortecimento_modais = (alpha_proporcional + omega_n_quadrado*beta_proporcional)/(2*omega_n);
fatores_amortecimento_modais = diag(fatores_amortecimento_modais);

% Construindo resposta em frequencia (H)
delta_frequencia = 0.1;
tempo_total = 1/delta_frequencia; % unidade em segundos
frequencia_amostragem = 44100/44; % 44100 amostras por segundo
frequencias_excitacao = [10:delta_frequencia:250]; % Hz
delta_tempo = 1/frequencia_amostragem;
tempos = [0:delta_tempo:tempo_total];
frequencias_excitacao = (0:length(tempos)-1)*frequencia_amostragem/length(tempos); 
omegas = 2*pi*frequencias_excitacao;
omegas_quadrado = omegas.^2;
H(numero_graus_liberdade, numero_graus_liberdade, length(omegas)) = 0;
for omega = 1:length(omegas)
	H(:,:,omega) = modos_vibracao* ... 
	(1./(omega_n_quadrado - omegas(omega).^2 + ... 
		i*2*diag(fatores_amortecimento_modais)*omega_n*omegas(omega)))*modos_vibracao';
end

H_sem_amortecimento(numero_graus_liberdade, numero_graus_liberdade, length(omegas)) = 0;
for omega = 1:length(omegas)
	H_sem_amortecimento(:,:,omega) = modos_vibracao* ... 
	(1./(omega_n_quadrado - omegas(omega).^2))*modos_vibracao';
end

%--
% Extraindo a resposta em 2 e excitação em 2
local_resposta = 2;
local_excitacao = 2;
H_22(1:length(omegas)) = 0;
H_22_sem_amortecimento(1:length(omegas)) = 0; 
for n = 1:length(omegas)
	H_22(n) = (H(local_resposta, local_excitacao, n));
	H_22_sem_amortecimento(n) = (H_sem_amortecimento(2,2,n)); 
end
% Calculando a acelerancia para resposta em 2 e excitação em 2
A_22 = omegas_quadrado.*H_22;
A_22(1) = mean(A_22);
A_22_sem_amortecimento = (omegas.^2).*H_22_sem_amortecimento;
A_22_sem_amortecimento(1) = mean(A_22_sem_amortecimento);

% Plotando Acelerancia numa Excitacao em 2 e Resposta Captada em 2
figure(1);
subplot(2,1,1);
semilogy(frequencias_excitacao, abs(A_22), 'black');
set(findobj(gca,'type','line'), 'LineWidth', 3);
hold on;
semilogy(frequencias_excitacao, abs(A_22_sem_amortecimento), 'blue');
title( ... 
'Modulo da Acelerancia numa Excitacao em 2 e Resposta em 2', ... 
'Interpreter','latex','FontSize',16);
xlabel('Frequencias [Hz]','Interpreter','latex','FontSize',16); 
ylabel('Magnitude','Interpreter','latex','FontSize',16);
legend('Com Amortecimento','Sem Amortecimento');
axis([10 250 min(abs(A_22_sem_amortecimento)) 1.1*max(abs(A_22_sem_amortecimento))]);
%
subplot(2,1,2);
plot(frequencias_excitacao, angle(A_22), 'black');
set(findobj(gca,'type','line'), 'LineWidth', 3);
hold on;
plot(frequencias_excitacao, angle(A_22_sem_amortecimento), 'blue');
title( ... 
'Fase da Acelerancia numa Excitacao em 2 e Resposta em 2', ... 
'Interpreter','latex','FontSize',16);
xlabel('Frequencias [Hz]','Interpreter','latex','FontSize',16); 
ylabel('Fase','Interpreter','latex','FontSize',16);
legend('Com Amortecimento','Sem Amortecimento');
axis([10 250 min(angle(A_22)) 1.1*max(angle(A_22_sem_amortecimento))]);

%--
% Extraindo a resposta em 4 e excitação em 2
local_resposta = 4;
local_excitacao = 2;
H_42(1:length(omegas)) = 0;
H_42_sem_amortecimento(1:length(omegas)) = 0;
for n = 1:length(omegas)
	H_42(n) = (H(local_resposta, local_excitacao, n));
	H_42_sem_amortecimento(n) = (H_sem_amortecimento(4,2,n));
end
% Calculando a acelerancia para resposta em 4 e excitação em 2
A_42 = omegas_quadrado.*H_42;
A_42(1) = mean(A_42);
A_42_sem_amortecimento = (omegas.^2).*H_42_sem_amortecimento;
A_42_sem_amortecimento(1) = mean(A_42_sem_amortecimento);
% Plotando Acelerancia numa Excitacao em 4 e Resposta Captada em 2
figure(2);
subplot(2,1,1);
semilogy(frequencias_excitacao, abs(A_42), 'black');
set(findobj(gca,'type','line'), 'LineWidth', 3);
hold on;
semilogy(frequencias_excitacao, abs(A_42_sem_amortecimento), 'blue');
title( ... 
'Modulo da Acelerancia numa Excitacao em 2 e Resposta em 4', ... 
'Interpreter','latex','FontSize',16);
xlabel('Frequencias [Hz]','Interpreter','latex','FontSize',16); 
ylabel('Magnitude','Interpreter','latex','FontSize',16);
legend('Com Amortecimento','Sem Amortecimento');
axis([10 250 min(abs(A_42_sem_amortecimento)) 1.1*max(abs(A_42_sem_amortecimento))]);
%
subplot(2,1,2);
plot(frequencias_excitacao, angle(A_42), 'black');
set(findobj(gca,'type','line'), 'LineWidth', 3);
hold on;
plot(frequencias_excitacao, angle(A_42_sem_amortecimento), 'blue');
title( ... 
'Fase da Acelerancia numa Excitacao em 2 e Resposta em 4', ... 
'Interpreter','latex','FontSize',16);
xlabel('Frequencias [Hz]','Interpreter','latex','FontSize',16); 
ylabel('Fase','Interpreter','latex','FontSize',16);
legend('Com Amortecimento','Sem Amortecimento');
axis([10 250 min(angle(A_42)) 1.1*max(angle(A_42_sem_amortecimento))]);

% ---------
% QUESTAO 5
frequencia_amostragem = 44100; % amostras por segundo
delta_tempo = 1/frequencia_amostragem;
frequencia_fundamental = 20; % unidades em Hz
periodo = 1/frequencia_fundamental;
tempos = [0:delta_tempo:periodo];
excitacao = tempos/(periodo/3);
menor_aproximado_1_3_periodo = abs(tempos - (periodo/3));
ponto_descontinuidade = find(menor_aproximado_1_3_periodo ...
== min(menor_aproximado_1_3_periodo));
excitacao(ponto_descontinuidade + 1:end) = 0;
excitacao_periodica = [excitacao excitacao excitacao];
tempos_total = linspace(0, 3*periodo, length(excitacao_periodica));
%figure(3);
%plot(tempos_total, excitacao_periodica);
% Montando a serie de fourier
x = linspace(0, 1, ponto_descontinuidade);
f = x;
N = 20; % numero de termos
fs = ones(1,length(x))/2;
for j = 1:N;
    fs_p = -1/pi*sin(2*j*pi*x)/j;
    fs = fs + fs_p;
end
excitacao_reconstruida = [fs excitacao(ponto_descontinuidade + 1:end)]; 
excitacao_reconstruida_periodica = [excitacao_reconstruida ...
 excitacao_reconstruida excitacao_reconstruida];

figure(6);
plot(tempos_total, excitacao_periodica,'LineWidth',3);
hold on
plot(tempos_total, excitacao_reconstruida_periodica, 'red','LineWidth',1);
title( ... 
'Excitacao de Deslocamento de Entrada no Sistema', ... 
'Interpreter','latex','FontSize',16);
xlabel('Tempo [s]','Interpreter','latex','FontSize',16); 
ylabel('Deslocamento','Interpreter','latex','FontSize',16);
legend('Sinal Original','Serie de Fourier do Sinal com 20 Termos');
hold off;
grid on;

