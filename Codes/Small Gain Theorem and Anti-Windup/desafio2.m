% ============================= Desafio 2 ==========================
clc;clear all; close all

%Modelos das plantas
L1 = 0.9; L2 = 0.7; L3 = 0.6; L4 = 0.4;
K1 = 1.3; K2 = 0.9; K3 = 1.2; K4 = .8;
tal1 = 1.2; tal2 = 1.1; tal3 = 0.8; tal4 = .9;

%Modelo Nominal
Ln = mean([L1 L2 L3 L4]);
Kn = mean([K1 K2 K3 K4]);
taln = mean([tal1 tal2 tal3 tal4]);

%==================> 1

w = logspace(-2,4,10000);
s = j*w;

%Modelo1
G1=K1.*exp(-L1*s)./(tal1*s+1);
Gn=Kn.*exp(-Ln*s)./(taln*s+1);
Delta1=abs(1-G1./Gn);
figure;
semilogx(w,Delta1);
title('Incertezas multiplicativas em função da frequência - Modelo 1','FontSize',16)
xlabel('Frequência w','FontSize',14);
ylabel('Incerteza Multiplicativa','FontSize',14)
grid on

%Modelo2
G2=K2.*exp(-L2*s)./(tal2*s+1);
Delta2=abs(1-G2./Gn);
figure;
semilogx(w,Delta2);
title('Incertezas multiplicativas em função da frequência - Modelo 2','FontSize',16)
xlabel('Frequência w','FontSize',14);
ylabel('Incerteza Multiplicativa','FontSize',14)
grid on

%Modelo3
G3=K3.*exp(-L3*s)./(tal3*s+1);
Delta3=abs(1-G3./Gn);
figure;
semilogx(w,Delta3);
title('Incertezas multiplicativas em função da frequência - Modelo 3','FontSize',16)
xlabel('Frequência w','FontSize',14);
ylabel('Incerteza Multiplicativa','FontSize',14)
grid on


%Modelo4
G4=K4.*exp(-L4*s)./(tal4*s+1);
Delta4=abs(1-G4./Gn);
figure;
semilogx(w,Delta4);
title('Incertezas multiplicativas em função da frequência - Modelo 4','FontSize',16)
xlabel('Frequência w','FontSize',14);
ylabel('Incerteza Multiplicativa','FontSize',14)
grid on

%==================> 2
close all
talc = 0.5;
Kc = (taln)/(Kn*(talc+Ln))
Ti = min(taln,4*(talc+Ln))

C = Kc*((s.*Ti+1))./(s.*Ti);


%==================> 3
Deltab = max([Delta1; Delta2; Delta3; Delta4],[],1);
Comp = (Gn.*C)./(1+Gn.*C);
mComp = abs(Comp);
figure;
semilogx(w,(Deltab.*mComp))
title('Função complementar de sensibilidade vezes limitante de incertezas multiplicativas','FontSize',16)
xlabel('Frequência w','FontSize',14);
ylabel('|C(jw)|Delta-Barra(w)','FontSize',14)
grid on











