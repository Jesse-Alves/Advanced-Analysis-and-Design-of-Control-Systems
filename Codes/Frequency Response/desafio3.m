% ============================= Desafio 3 ==========================
clc;clear all; close all

%w = logspace(-2,4,10000);
%==================> 1
s = 0;
G0 = (0.5)/((s^2+0.6*s+1)*(0.1*s+1));
Kb = 1/G0

%==================> 2
%Planta
s = tf([1 0],1);
G = (0.5)/((s^2+0.6*s+1)*(0.1*s+1));

P = Kb*G;
Largura_Banda = bandwidth(P)

%%
%==================> 3
Kc = 1;
G1 = Kc*P;
[Gm1,Pm1,Wcg1,Wcp1] = margin(G1);
display('A margem de Fase do Sisitema G1:')
Pm1

display('Valor do Fi máximo:')
fi_m = 60 - Pm1 + 12

syms alfa
a = solve(sind(fi_m) == (1-alfa)/(1+alfa));
display('Valor de Alfa:')
a = double(a)

sysComp = G1*(1/sqrt(a));
[Gm,Pm,Wcg,Wcp] = margin(sysComp);

display('A frequencia de cruzamento de ganho:')
Wcp
display('Valor de T:')
T = 1/(sqrt(a)*Wcp)
display('O controlador:')
C = (Kc*(T*s+1))/(a*T*s+1)

FTma = C*P;
figure
margin(FTma)
[Gm,Pm,Wcg,Wcp] = margin(FTma);
display('Margem de Fase da Malha Aberta:')
Pm
grid on


%%
%==================> 4
Kc = 1;
G1 = Kc*P;
[Gm1,Pm1,Wcg1,Wcp1] = margin(G1);
display('A margem de Fase do Sisitema G1:')
Pm1

display('Valor do Fi máximo:')
fi_m = 60 - Pm1 + 24 %Folga de 24 graus

syms alfa
a = solve(sind(fi_m) == (1-alfa)/(1+alfa));
display('Valor de Alfa:')
a = double(a)

sysComp = G1*(1/sqrt(a));
[Gm,Pm,Wcg,Wcp] = margin(sysComp);

display('A frequencia de cruzamento de ganho:')
Wcp
display('Valor de T da Questão 4:')
T = 1/(sqrt(a)*Wcp)
display('O controlador da Questão 4:')
C = (Kc*(T*s+1))/(a*T*s+1)

FTma = C*P;
figure
margin(FTma)
[Gm,Pm,Wcg,Wcp] = margin(FTma);
display('Margem de Fase da Malha Aberta:')
Pm
grid on

%%
%==================> 5
%Identificar o valor do ganho na frequencia w = 2,5
bode(FTma)
grid on

w = 2.5;
s = j*w;
syms kc
FTma = ((kc*(T*s+1))/(a*T*s+1))*Kb*((0.5)/((s^2+0.6*s+1)*(0.1*s+1)));
Kc = solve(abs(FTma) == 1/sqrt(2));
display('Kc para largura de banda de 2,5 1/rad*s')
Kc = double(Kc) 

%Sistema com Ganho Kc modificado
s = tf([1 0],1);
display('O controlador:')
C = (Kc*(T*s+1))/(a*T*s+1)
FTma = C*P;
[Gm,Pm,Wcg,Wcp] = margin(FTma);
display('A Nova Margin de fase para Kc = 1.4297:')
Pm

%%
%==================> 6
Kc = 1.4297; %<===================     NOVO K
G1 = Kc*P;
[Gm1,Pm1,Wcg1,Wcp1] = margin(G1);
display('A margem de Fase do Sisitema G1:')
Pm1

display('Valor do Fi máximo:')
fi_m = 60 - Pm1 + 24 %Folga de 24 graus

syms alfa
a = solve(sind(fi_m) == (1-alfa)/(1+alfa));
display('Valor de Alfa:')
a = double(a)

sysComp = G1*(1/sqrt(a));
[Gm,Pm,Wcg,Wcp] = margin(sysComp);

display('A frequencia de cruzamento de ganho:')
Wcp
display('Valor de T da Questão 6:')
T = 1/(sqrt(a)*Wcp)
display('O controlador da Questão 6:')
C = (Kc*(T*s+1))/(a*T*s+1)

FTma = C*P;
figure
margin(FTma)
[Gm,Pm,Wcg,Wcp] = margin(FTma);
display('Margem de Fase da Malha Aberta:')
Pm
display('E a Banda de Passagem:')
bandwidth(FTma)
grid on

%%
%==================> 8
% Para C4
s = 0;
Kc = 1;
T4 = 1.282;
Tat = 5*T4;
C4= (1.282*s + 1)/(0.2698*s + 1);

syms beta
b = solve(C4*Kb*((beta*(Tat*s+1))/(beta*Tat*s+1))*G0 == 5);
b = double(b)

s = tf([1 0],1);
C4= (1.282*s + 1)/(0.2698*s + 1);
display('Controlador do item 4 Modificado:')
C4new = C4*Kb*((b*(Tat*s+1))/(b*Tat*s+1))

%%
%==================> 8
% Para C6 
s = 0;
Kc = 1.4297;
T6 = 1.3711;
Tat = 5*T6;
C6 = (1.96*s + 1.43)/(0.143*s + 1);

syms beta
b = solve(C6*Kb*((beta*(Tat*s+1))/(beta*Tat*s+1))*G0 == 5);
b = double(b)

s = tf([1 0],1);
C6 = (1.96*s + 1.43)/(0.143*s + 1);
display('Controlador do item 6 Modificado:')
C6new = C6*Kb*((b*(Tat*s+1))/(b*Tat*s+1))


