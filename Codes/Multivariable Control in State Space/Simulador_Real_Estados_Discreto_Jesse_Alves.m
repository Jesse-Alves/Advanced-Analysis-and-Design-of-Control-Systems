clc
%clear all
close all

%% Define caso (0 - nominal, 1 -caso incerto)
flag=0;

% Variável complexa s (função de transferência)
s=tf([1 0],1);

%% Sistema simulado
P = G;            % Processo contínuo
L=0.0; % atraso

%% Modelo nominal para fins de projeto
Pn=P;
Ln=0; % atraso


%% Tempo total de simulação 
tsim=100; % Tempo de simulação

%% Tempo de alteração da
tset=0;  % Instante de alteração da referência

%% Define tipo de perturbação (0 - degrau, 1-rampa, 2 senoide)

tipo_u=0;   % tipo de perturbação da entrada
tipo_y=0;   % tipo de perturbação da saída

tpu=40;     % Intante inicial da perturbação da entrada
du=80;      % Intervalo de duração da perturbação da entrada
Apu=0.1;   % Amplitude da perturbação da entrada

tpy=20;     % Intante inicial da perturbação da saída
dy=80;       % Duração da perturbação da saída
Apy=0.2;    % Amplitude da perturbação da entrada

omega_u=2;  % Frequêcia da perturação da entrada (apenas para tipo 2)
omega_y=2;  % Frequêcia da perturação da saída   (apenas para o tipo 2)

%% Controlador 

Ts=ta;                   % Período de amostragem
Pnd= ssGd;                     % Modelo discretizado

tipo_cont = 1;
if tipo_cont==0 
    % ==================> SEM INTEGRADOR
    Ad=Pnd.A;
    Bd=Pnd.B;
    Cd=Pnd.C;
    
    Ka = Kd;

else
    % ==================> COM INTEGRADOR
    Ad=Pnd.A;
    Bd=Pnd.B;
    Cd=Pnd.C;    
    
    no=size(Cd,1);
    ni=size(Bd,2);
    ns=size(Ad,2);    
    
    Aa=[Ad zeros(ns,no);Cd*Ad eye(no)];
    Ba=[Bd; Cd*Bd];
    Ca=[Cd zeros(no,no); zeros(no,ns) eye(no)];

    Ka= Ka;

    % Elemento para Estimação de estados
    %La=place(Aa',Ca',0.4*ones(size(Aa,1),1))'; %%%%%%%%

    Ki = Ki;
    
end
%% Define passo de integração e passo de armazenamonto 
dT=Ts;     % passo de amostragem 
h=Ts/100;  % passo de integração

%%----------------- finaliza definição de parâmetros ------------------------

if flag==0
  P=Pn;
  L=Ln;
end

% Representação em espaço de estados (acelera a simulação)
ssP= ssG; Ap=ssP.a; Bp=ssP.b; Cp=ssP.c; Dp=ssP.d;
ssPn=ss(Pn); Apn=ssPn.a; Bpn=ssPn.b; Cpn=ssPn.c; Dpn=ssPn.d;


% Arredonda o passo de integração tal que dT seja múltiplo de h
m=ceil(dT/h);
h=dT/m;

% Atraso de tempo discretro
Ld=round(L/h);
U_atraso=zeros(1,Ld+1);


display('Atraso múltiplo do passo - atraso considerado:')
display(Ld*h)

% Calcula número de passos de simulação
kfinal=ceil(tsim/dT);


% Estabelece condições iniciais dos elementos dinâmicos
xp0=zeros(size(Ap,1),1);
xhat=zeros(size(Aa,1),1);
up=[0;0];

% Simula sistema de controle
 
for k=1:kfinal
  
   t(k)=(k-1)*dT;  % marca o tempo  

   % Estabelece as condições iniciais (apenas no primeiro passo)
   
   if k==1
    
    xp=xp0;   % condição inicial do processo
    y=Cp*xp0;  % condição inicial da saída  (assumindo Dp=0)  
    yp=y;
 
   end
 

 % Define mudança de referência
  if k<tset/dT+1
    r(:,k)=[0;0];
  else
    r(:,k)=[2;1]; 
  end
   

 % Define perturbação na entrada do processo
  if k<=tpu/dT-1   
      qu=0;
  elseif k>tpu/dT&& k<= (tpu+du)/dT -1
      
         if  tipo_u==0           
                  qu=Apu;
         elseif tipo_u==1
                  qu=qu+Apu*dT;
         else
                  qu=sin(omega_u*k*dT);
         end
  else
         if tipo_u~=1 
             qu=0;
          end 
  end
        
 % Define perturbação na saída do processo
  if k<=tpy/dT-1   
      qy=0;
  elseif k>tpy/dT&& k<= (tpy+dy)/dT -1
      
         if  tipo_y==0           
                  qy=Apy;
         elseif tipo_y==1
                  qy=qy+Apy*dT;
         else
                  qy=sin(omega_y*k*dT);
         end
  else           
          if tipo_y~=1 
             qy=0;
          end 
  end    
  
%   % SEM INTEGRADOR
%   delta_u = -Ka*xhat;
%   u=delta_u+up;
%   up=u;
%   delta_y = y-yp;
%   y_e=xhat(end-no+1:end);
%   delta_ye=Cd*xhat(1:end-no);
%   
%   xhat=Aa*xhat+La*[delta_y-delta_ye;y-y_e]+Ba*delta_u;
%   
%   yp=y;
  
  
  % COM INTEGRADOR
  delta_u=-Ka*xhat+Ki*r(:,k);
  u=delta_u+up;
  up=u;
  delta_y=y-yp;
  %xhat=(Aa-La*Ca)*xhat+La*[delta_y;y]+Ba*delta_u;
  y_e=xhat(end-no+1:end);
  delta_ye=Cd*xhat(1:end-no);
% xhat=Aa*xhat+La*[delta_y-delta_ye;y-y_e]+Ba*delta_u; %Com estimador
  xhat=Aa*xhat + Ba*delta_u;
  yp=y;
   
  
   for l=1:m      

     %Calcula saída do processo
     U_atraso=[u U_atraso(:,1:Ld)];
     u_atraso=U_atraso(:,Ld+1);
     
   % [y_aux,xp]=sim_euler(Ap,Bp,Cp,Dp,xp,h,uatraso+qu);   
     xp=h*Ap*xp+xp+h*Bp*(u_atraso+qu);
     y_aux=Cp*xp+Dp*(u_atraso+qu);
     y=y_aux+qy;
     
    
     
     % Armazena valores calculados
     if l==1
       
       Xp(:,k)=xp;   
         
       U(:,k)=u_atraso;
       Y(:,k)=y;
       Qy(:,k)=[qy;qy];
       Qu(:,k)=[qu;qu];
    
    end
              
   end
       
end


nfont=15;   % define tamanho da fonte de texto
nlinha=3;   % define espessura da linha

figure(1)

%t = linspace(0,tsim,length(r(1,:)));
subplot(2,1,1)
plot(t,r(1,:),'--r','linewidth',nlinha-0.5)
hold on
plot(t,Y(1,:),'k','linewidth',nlinha)
plot(t,r(2,:),'--r','linewidth',nlinha-0.5)
plot(t,Y(2,:),'b','linewidth',nlinha)
grid on
axis([0 t(end) min(min(Y))-0.2  max(max(Y))+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Saída');
%hl=legend('r(t)', 'y(t)');
%set(hl,'fontsize',nfont)
%set(hl,'location','southeast')
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

subplot(2,1,2)
plot(t,U(1,:),'k','linewidth',nlinha)
hold on
plot(t,U(2,:),'b','linewidth',nlinha)
grid on
axis([0 t(end) min(min(U))-0.2 max(max(U))+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Controle');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)


figure(2)

subplot(2,1,1)
plot(t,Qy(1,:),'b','linewidth',nlinha)
hold on
plot(t,Qy(2,:),'--r','linewidth',nlinha)
grid on
axis([0 t(end) min(min(Qy))-0.2  max(max(Qy))+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Perturbação - Saída');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

subplot(2,1,2)
plot(t,Qu(1,:),'b','linewidth',nlinha)
hold on
plot(t,Qu(2,:),'--r','linewidth',nlinha)
grid on
axis([0 t(end) min(min(Qu))-0.2 max(max(Qu))+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Perturbação - Entrada');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)
