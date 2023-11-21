%clc
%clear all
% ================ SIMULA��O DISCRETO ================
%tempo de simulação
tsim=100;

% sistema em malha aberta
G = tf([-0.2 2],[1 2 1]);
Gs=ss(G);
Ag=Gs.A; Bg=Gs.B; Cg=Gs.C; Dg=Gs.D;

% Período de amostragem
Ts=0.2;
% Razão-> Período de amostragem/passo de integração
M=20;

h=Ts/M;

% controlador

C=tf([b1 b0],[a1 a0]);
Cs=ss(C);
Ac=Cs.A; Bc=Cs.B; Cc=Cs.C; Dc=Cs.D;

Cd=c2d(C,Ts,'tustin');

[numC,denC]=tfdata(Cd,'v');

% Normaliza den(1) - não é necessário (fim didático)
numC=numC/denC(1);
denC=denC/denC(1);

% número de passos (loop de integração)

Nsim=tsim/h;


% referência
tref=10;
Nref=tref/h;
ref=[zeros(Nref,1); ones(Nsim-Nref,1)];


% perturbação controle
tqu=40;
Nqu=tqu/h;
qu=[zeros(Nqu,1); 0.2*ones(Nsim-Nqu,1)];

% perturbação saída
tqy=70;
Nqy=tqy/h;
qy=[zeros(Nqy,1); 0.2*ones(Nsim-Nqy,1)];

% consições iniciais

xg=zeros(size(Ag,1),1);
yg=Cg*xg;
y=yg;

xc=zeros(size(Ac,1),1);

ep=0;
up=0;

for k=1:Nsim
  
 % Neste código, u(kTs) é atualiadado a cada a cada múltiplo M de h.
 % Em outras palavras, Ts=M*h sendo h o passo de integração e Ts o período de amostragem. 
  
  if mod(k,M)==1 % Atualiza quanto resto da divisão de k por M igual a 1 => atualize u(kTs)
  
   e=ref(k)-y;
   
  % Dinâmica Controlador PI de Tempo Discreto
  
  u=-denC(2)*up+numC(1)*e+numC(2)*ep;
  
  % Armazena o controle passado (u(kTs)) e o erro passado (e(kTs)) para uso em (k+1)Ts
  up=u;
  ep=e;
 
  end
  
  % Dinâmica Processo
    
  xg=xg+h*Ag*xg+h*Bg*(u+qu(k));
  yg=Cg*xg+Dg*u;
  
  y=yg+qy(k);
  
  % "Scope"
  
  tempo(k)=k*h;
  
  R(k)=ref(k);
  Y(k)=y;
  U(k)=u;
  
  
  
end

figure(1)
plot(tempo,R,'--r','LineWidth',3)
hold on
plot(tempo,Y,'k','LineWidth',3)
grid on

% figure(2)
% plot(tempo,U,'b','LineWidth',3)
% grid on









