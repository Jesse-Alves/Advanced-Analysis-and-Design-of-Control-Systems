%% ============================= Desafio 7 ==========================
clc;clear all; close all

%pkg load control

%% Questão 1 e 2

disp('Sistema Linearizado')
syms x1 x2
A = [0 1; -1-0.6*x1*x2 -0.3*x1^2+0.3]
B = [0;1]

disp('Modelo 1')
x1 = 1;
x2 = 0;
A1 = [0 1; -1-0.6*x1*x2 -0.3*x1^2+0.3]
B1 = [0;1]

disp('Modelo 2')
x1 = 4;
x2 = 0;
A2 = [0 1; -1-0.6*x1*x2 -0.3*x1^2+0.3]
B2 = [0;1]

%% Questão 3
auto_val1 = eig(A1)
auto_val2 =eig(A2)

%% Questão 4
R = 0.1;
Q = eye(2);

disp('Ganho K para o Modelo 1')
[K1,S,E] = lqr(A1, B1, Q, R);
K1

disp('Ganho K para o Modelo 1')
[K2,S,E] = lqr(A2, B2, Q, R);
K2

%% SIMULAÇÃO TEMPORAL
dT = 0.01;
tsim = 10;
t = 0:dT:tsim;

kfinal = length(t);

%Valores iniciais
x0 = [-5;-5]; %Cond. Inicial
x = zeros(2,kfinal);
x(:,1) = x0;
u = zeros(1,kfinal);

%===================== LEIS DE CONTROLE =========================
Caso = 1;
if Caso == 0
    K = [0 0];
    ueq = 0; 
    xeq = [0; 0];
elseif Caso == 1
    K = K1;
    ueq = 1; 
    xeq = [1; 0];
elseif Caso == 2
    K = K2;
    ueq = 4; 
    xeq = [4; 0];
end

% TENTAR INVERTER OS Ks

%================================================================
%Laço de Simulação
for k = 1:kfinal-1
    u(k) = -K*(x(:,k) - xeq) + ueq;
    dx = dT*[x(2,k); -x(1,k)+0.3*(1-x(1,k)^2)*x(2,k) + u(k)];
    x(:,k+1) = x(:,k) + dx;
end

%Plotagem
figure;
plot(t,x(1,:),'k','linewidth',2)
hold on
plot(t,x(2,:),'b','linewidth',2)
grid on

%Linha de ref para o estado 1
plot(t,xeq(1)*ones(1,kfinal),'r--','linewidth',2)
%Linha de ref para o estado 1
plot(t,zeros(1,kfinal),'r--','linewidth',2)

xlabel('Tempo (s)','FontSize',16)
ylabel('Estados do Sistema - Oscilador de Van der Pol','FontSize',16)
legend({'x_1(t)','x_2(t)'}, 'FontSize',15)
