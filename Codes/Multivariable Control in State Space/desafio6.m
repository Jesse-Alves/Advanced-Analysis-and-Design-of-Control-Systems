%% ============================= Desafio 6 ==========================
clc;clear all; close all

%pkg load control

s = tf([1 0],1);
%Planta
g1 = 2/(10*s+1);
g2 = 0.8/((10*s+1)*(2*s+1));
g3 = 0.6/((10*s+1)*(2*s+1));
g4 = 1.5/(10*s+1);
G = [g1 g2;g3 g4];


% ===> Perido de Amostragem
[y,tempo] = step(g1);
vet = find(y >= 0.98*y(length(y)));
tempo(vet(1))
ta1 = (tempo(vet(1)))/20;

[y,tempo] = step(g2);
vet = find(y >= 0.98*y(length(y)));
tempo(vet(1))
ta2 = (tempo(vet(1)))/20;

[y,tempo] = step(g3);
vet = find(y >= 0.98*y(length(y)));
ta3 = (tempo(vet(1)))/20;

[y,tempo] = step(g4);
vet = find(y >= 0.98*y(length(y)));
ta4 = (tempo(vet(1)))/20;

disp('O menor valor que atende a todas as FTs')
t_min = max([ta1 ta2 ta3 ta4])

disp('Arbitra-se portanto um tempo de amostragem')
ta = 2.07

%% Controlador u_eq
yeq = [2;1];

s = 0;
g1 = 2/(10*s+1);
g2 = 0.8/((10*s+1)*(2*s+1));
g3 = 0.6/((10*s+1)*(2*s+1));
g4 = 1.5/(10*s+1);
G0 = [g1 g2;g3 g4]; 

disp('Vetor Controle de Equilibrio')
ueq = inv(G0)*yeq

%% Espaço de Estados
s = tf([1 0],1);
ssGa = ss(G);
Aa = ssGa.A;
Ba = ssGa.B;
Ca = ssGa.C;

ssG = minreal(ssGa);
disp(' ')
disp('Matrizes no Espaço de Estados no em Tempo Continuo')
disp(' ')
A = ssG.A
B = ssG.B
C = ssG.C

disp(' ')
disp('G(s) = C(sI - A)^-1B')
disp(' ')
G_calc = Ca*inv(s*eye(length(Aa)) - Aa)*Ba

%% Estado de Equilibrio - Tempo Contínuo
disp('Vetor Estado de Equilibrio - Tempo Continuo')
xeq = -inv(A)*B*ueq

%% Estado de Equilibrio Discreto
z = tf([1 0],1,ta);

%Sistema sem realizacao minima
ssGad = c2d(ssGa,ta);
Aad = ssGad.A;
Bad = ssGad.B;
Cad = ssGad.C;

ssGd = c2d(ssG,ta);
disp(' ')
disp('Matrizes no Espaço de Estados em Tempo Discreto')
disp(' ')
Ad = ssGd.A
Bd = ssGd.B
Cd = ssGd.C
Dd = ssGd.D

disp(' ')
disp('Gd(z) = Cd(zI - Ad)^-1Bd')
disp(' ')
Gd_calc = Cad*inv(z*eye(length(Aad)) - Aad)*Bad

%% Estado de Equilibrio - Tempo Discreto
disp('Vetor Estado de Equilibrio - Tempo Discreto')
xeq_d = inv(eye(length(Ad)) - Ad)*Bd*ueq

disp(' ')
disp('Verificando igualdade xeq=xeq_d')
disp(' ')
xeq
xeq_d

disp(' ')
disp('yeq igual ao pedido')
disp(' ')
yeq_calc = Cd*xeq_d

%% Calculo do ganho K

lambda = exp(-ta/4);    
polos = [lambda lambda lambda lambda];
%polos = [lambda lambda lambda-0.1 lambda-0.1];
disp('Ganho Kd - Sem Integrador')
Kd = place(Ad,Bd,polos)

%% Calculo do ganho K - COM INTEGRADOR
Aa = [Ad zeros([size(Ad,1),size(Bd,2)]);-Cd zeros(size(Cd,1))];
Ba = [Bd;-Dd];

polos_a = [lambda lambda lambda lambda lambda lambda];
%Ganho K Aumentado
Ka = place(Aa,Ba,polos_a);

[n,m] = size(Bd);
[q,n] = size(Cd);
disp('Ganho de Realimentação de Estado')
Kd_i = Ka(:,1:n)
disp('Ganho do Integrador')
Ki = Ka(:,n+1:n+q)

%% CHAMAR SIMULAÇÃO
%Simulador_Real_Estados_Discreto_Jesse_Alves()
%Simula_MIMO()
