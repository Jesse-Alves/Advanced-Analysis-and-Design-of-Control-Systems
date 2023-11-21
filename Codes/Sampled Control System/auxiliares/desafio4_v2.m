%% ============================= Desafio 4 ==========================
clc;clear all; close all

%Planta
s = tf([1 0],1);
Gs = (0.2*(10-s))/((s+1)^2);

%Tempo de Acomodação Malha Aberta 
inform_5 = stepinfo(Gs,'SettlingTimeThreshold',0.05);
inform_2 = stepinfo(Gs);

ts_5 = inform_5.SettlingTime;
ts_2 = inform_2.SettlingTime;

disp('Período de Amostragem')
ta = (ts_5)/20

disp('Tempo de Acomodação Desejado')
ts = (ts_2)/2 

%% ALOCAÇÃO DE POLOS
z = tf([1 0],1,ta);

% ===> Criterios de Projeto
%Polos repetidos => ts para 2% = 6*tau
tau = ts/6;

% ===> Polos do sistema no contínuo
s1 = -1/tau;
s2 = -1/tau;
s3 = -30; % Outros dois polos 15x maiores
s4 = -30;

% ===> Polos do sistema no discreto
z1 = exp(ta*s1);
z2 = exp(ta*s2);
z3 = exp(ta*s3); 
z4 = exp(ta*s4);

%Polinomio Caracteristico Desejado
D0 = poly([z1 z2 z3 z4]);

%Planta Discretizada
Gz = c2d(Gs,ta,'zoh');
[numGz denGz] = tfdata(Gz,'v');

%Comparando polinomios na forma matricial

%Fazendo A0 = 0 por conta do número de equações e variáveis
A0 = 0;
M = [numGz(3) 0 0 0 0;
     numGz(2) denGz(3) numGz(3) 0 0;
     0 denGz(2) numGz(2) denGz(3) numGz(3);
     0 1 0 denGz(2) numGz(2);
     0 0 0 1 0];
 
coef_controlador = inv(M)*(fliplr(D0))';
B0 = coef_controlador(1);
A1 = coef_controlador(2);
B1 = coef_controlador(3);
A2 = coef_controlador(4);
B2 = coef_controlador(5);

disp('Controlador por Alocação de Polos')
Cz = (B0 + B1*z + B2*z^2)/(A0 + A1*z + A2*z^2)
[numCz denCz] = tfdata(Cz,'v');

disp('Filtro de Referência')
Fz = tf(poly([z3 z4]),numCz,ta);
Hz = feedback(Cz*Gz,1);
K = evalfr(Fz*Hz,1);

K = 1/K;
Fz = K*Fz;

[numFz denFz] = tfdata(Fz,'v');
%step(Fz*Hz)

%Teste_Simula_Discreto()
%% CONTROLADOR IMC
%Fr_til = 0.4^2/(z - 0.6)^2;
%B_mais = 0.011856;
%B_menos = z - 

%Cz_IMC = Hd_til*Az/(1-Hd_til*B_z);
%[numIMC denIMC] = tfdata(Cz_IMC,'v');

C_IMC = (0.02057/0.011856)*(z-0.7853)^2/((z-0.6)^2-0.02057*(z+6.779));
[numC_IMC denC_IMC] = tfdata(C_IMC,'v');












