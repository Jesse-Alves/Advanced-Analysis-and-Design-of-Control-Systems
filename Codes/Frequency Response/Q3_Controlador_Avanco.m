%Especificações 
clc;clear all;
MF = 60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Função de Transferência:
s = tf('s');
K = 2;
%1.4298
G = 0.5/((s^2 + 0.6*s + 1)*(0.1*s + 1));
[num , den] = tfdata(G, 'v');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Encontrando a Margem de Fase de KG(s):
sys = double(K)*G;
[Gm,Pm,Wcg,Wcp] = margin(sys);

fprintf('\nPara KG(s) tem-se: MF = %f e wcg = %f; MG = %f e wcf = %f\n',Pm,Wcp,Gm,Wcg);

%Determinando alfa a partir de phi_m:
ajuste_fase = 12;
phi_m = MF - Pm + ajuste_fase;
alfa = (1 + sind(phi_m))/(1 - sind(phi_m));

fprintf('\nPhi_m = %f e alfa = %f\n',phi_m,alfa);

%Determinando T a partir de wcg do sistema compensado:
sys2 = sqrt(alfa)*sys;
[Gm2,Pm2,Wcg2,Wcp2] = margin(sys2);

T = 1/(Wcp2*sqrt(alfa));

fprintf('\nFrequência de cruzamento de ganho do sistema compensado: Wcg = %f -> T = %f e alfa*T = %f\n',Wcp2,T,alfa*T);

%Determinando o controlador em Avanço e fazendo a verificação do sistema:
CAV_s = double(K)* (alfa*T*s + 1)/(T*s + 1);

sys3 = CAV_s*G;
[Gm3,Pm3,Wcg3,Wcp3] = margin(sys3);
wb = bandwidth(sys3);
malha_fechada = feedback(sys3,1);
info = stepinfo(malha_fechada);

fprintf('\nO controlador sintetizado é:\n');
CAV_s
fprintf('\nC(s)G(s) = ');
sys3
fprintf('A margem de ganho, fase, frequência de cruzamento de ganho do sistema compensado e a largura de banda são, respectivamente: MG = %f, MF = %f , wcg = %f, wb = %f\n',Gm3,Pm3,Wcp3,wb);
fprintf('\nA função de transferência de malha fechada é:\n');
malha_fechada
fprintf('O máximo sobressinal e o tempo de acomodação são, respectivamente: Mp = %f e ts = %f\n', info.Overshoot, info.SettlingTime);
