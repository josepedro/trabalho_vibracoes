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
c_2 = 120; % unidade em N.s/m
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
c_4 = alpha_proporcional*M(1,1) + beta_proporcional*K(1,1) - c_1 - c_2;
c_5 = alpha_proporcional*M(3,4) + beta_proporcional*K(3,4);
c_6 = alpha_proporcional*M(2,2) + beta_proporcional*K(2,2) - c_3 - c_2;
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

% Construindo a forca de excitacao
delta_frequencia = 0.1;
tempo_total_excitacao = 1/delta_frequencia; % unidade em segundos
frequencia_amostragem = 44100/44; % 44100 amostras por segundo
frequencias_excitacao = [10:delta_frequencia:250]; % Hz
delta_tempo = 1/frequencia_amostragem;
tempos = [0:delta_tempo:tempo_total_excitacao];
forca_excitacao = 0;
% Realizando o somatorio das forcas
for frequencia = 1:length(frequencias_excitacao)
	forca_excitacao = forca_excitacao + sin(frequencias_excitacao(frequencia)*2*pi*tempos);
end
% Plotando a banda de frequencia de excitacao da forca
frequencias_ = (0:length(forca_excitacao)-1)*frequencia_amostragem/length(forca_excitacao); 
figure(1);
plot(frequencias_, abs(fft(forca_excitacao)));
title('Espectro de Frequencias da Forca de Excitacao','Interpreter','latex','FontSize',16);
xlabel('Frquencias [Hz]','Interpreter','latex','FontSize',16); 
ylabel('Amplitude','Interpreter','latex','FontSize',16);
axis([0 (frequencia_amostragem/2) min(abs(fft(forca_excitacao))) ... 
	max(abs(fft(forca_excitacao)) + 0.1*abs(fft(forca_excitacao)))]);


% Construindo resposta em frequencia (H)
omegas = 2*pi*frequencias_;
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

% Extraindo a resposta em 2 e excitação em 2
local_aplicacao_forca = 2;
local_captacao_resposta = 2;
H_22(1:length(omegas)) = 0;
H_22_sem_amortecimento(1:length(omegas)) = 0; 
for n = 1:length(omegas)
	H_22(n) = abs(H(local_captacao_resposta,local_aplicacao_forca,n));
	H_22_sem_amortecimento(n) = abs(H_sem_amortecimento(2,2,n)); 
end

figure;
semilogy(frequencias_, H_22, 'black');
set(findobj(gca,'type','line'), 'LineWidth', 3);
hold on;
semilogy(frequencias_, H_22_sem_amortecimento, 'blue');
axis([0 250 0 0.0001])


% Extraindo a resposta em 4 e excitação em 2
H_42(1:length(omegas)) = 0;
H_42_sem_amortecimento(1:length(omegas)) = 0;
for n = 1:length(omegas)
	H_42(n) = abs(H(4,local_aplicacao_forca,n));
	H_42_sem_amortecimento(n) = abs(H_sem_amortecimento(4,2,n));
end

figure;
semilogy(frequencias_, H_42, 'black');
set(findobj(gca,'type','line'), 'LineWidth', 3);
hold on;
semilogy(frequencias_, H_42_sem_amortecimento, 'blue');
axis([0 250 0 0.0001])

