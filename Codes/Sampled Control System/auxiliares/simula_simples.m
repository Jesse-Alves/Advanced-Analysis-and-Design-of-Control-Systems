clc
clear 

%tempo de simulação
tsim=100;

%passo de simulação
h=0.01;


% sistema em malha aberta
G=tf(1,[1 1])*tf(1,[0.1 1]);
Gs=ss(G);
Ag=Gs.A; Bg=Gs.B; Cg=Gs.C; Dg=Gs.D;

% controlador

C=tf(1*[1 1],[1 0]);
Cs=ss(C);
Ac=Cs.A; Bc=Cs.B; Cc=Cs.C; Dc=Cs.D;

% número de passos (loops)

Nsim=tsim/h;


% referência

tref=10;
Nref=tref/h;
ref=[zeros(Nref,1); ones(Nsim-Nref,1)];


% perturbação controle
tqu=40;
Nqu=tqu/h;
qu=[zeros(Nqu,1); 0.3*ones(Nsim-Nqu,1)];

% perturbação saída
tqy=70;
Nqy=tqy/h;
qy=[zeros(Nqy,1); -0.3*ones(Nsim-Nqy,1)];

% consições iniciais

xg=zeros(size(Ag,1),1);

yg=Cg*xg;
y=yg;

xc=zeros(size(Ac,1),1);


for k=1:Nsim
  
  e=ref(k)-y;
  
  
  % Dinâmica Controlador
  
  xc=xc+h*Ac*xc+h*Bc*e;
  yc=Cc*xc+Dc*e;
  
  u=yc;
  
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


% gera gráficos

figure(1)
plot(tempo,R,'linewidth',3,'-r')
hold on
plot(tempo,Y,'linewidth',3,'--b')

figure(2)
plot(tempo,U,'linewidth',3,'-r')
