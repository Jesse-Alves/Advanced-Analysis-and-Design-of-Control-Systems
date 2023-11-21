%% ============================= Desafio 4 ==========================
clc;clear all; close all

%Planta
s = tf([1 0],1);
Gs = (0.2*(10-s))/((s+1)^2);

%Tempo de acomodação de 2%
[y,tempo] = step(Gs);
vet = find(y >= 0.98*y(length(y)));
ts_2 = tempo(vet(1));


%Tempo de acomodação de 5%
vet = find(y >= 0.95*y(length(y)));
ts_5 = tempo(vet(1));

disp('Período de Amostragem')
ta = (ts_5)/20

disp('Tempo de Acomodação Desejado')
ts = 3 

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
D0 = [0 D0];

%Planta Discretizada
Gz = c2d(Gs,ta,'zoh');
[numGz denGz] = tfdata(Gz,'v');

%Comparando polinomios na forma matricial

%Fazendo A0 = 0 por conta do número de equações e variáveis
%A0 = 0;
M = [denGz(3) numGz(3) 0 0 0 0;
     denGz(2) numGz(2) denGz(3) numGz(3) 0 0;
     1 0 denGz(2) numGz(2) denGz(3) numGz(3);
     0 0 1 0 denGz(2) numGz(2);
     0 0 0 0 1 0;
     -1 0 -1 0 -1 0];
 
coef_controlador = inv(M)*(fliplr(D0))';
A0 = coef_controlador(1);
B0 = coef_controlador(2);
A1 = coef_controlador(3);
B1 = coef_controlador(4);
A2 = coef_controlador(5);
B2 = coef_controlador(6);

disp('Controlador por Alocação de Polos')
Cz = (B0 + B1*z + B2*z^2)/(A0 + A1*z + A2*z^2)
[numCz denCz] = tfdata(Cz,'v');

disp('Filtro de Referência')
Fz = tf(poly([z3 z4]),numCz,ta);
%Ajustando o ganho
K = evalfr(Fz*feedback(Cz*Gz,1),1);
Fz = (1/K)*Fz
[numFz denFz] = tfdata(Fz,'v');

disp('Indices de desempenho do Controlador por Alocacao de Polos')
stepinfo(Fz*feedback(Cz*Gz,1))
%step(Fz*feedback(Cz*Gz,1))
%grid on

%% CONTROLADOR IMC

%Planta Discretizada 
Pn = Gz;
B_ = z + 6.779; %Parcela que não pode ser cancelada

Fr_til = (0.4^2)/((z-0.6)^2); %Arbitrada para ts ser 3s 
[numFr_til denFr_til] = tfdata(Fr_til,'v');
%stepinfo(Fr_til)

disp('Filtro de Robustez')
Fr = B_*Fr_til;
Fr = (1/dcgain(Fr))*Fr_til*B_
[numFr denFr] = tfdata(Fr,'v');
%dcgain(Fr) %Ganho unitario

%%
disp('FT do Controlador IMC')
C_IMC = (1.7348*(z-0.7853)^2)/((z-1)*(z-0.22));

disp('Indices de desempenho do Controlador IMC')
stepinfo(feedback(C_IMC*Pn,1))


%Chamar Simulação
Simulador_Controle_Discreto_Jesse_D4()











