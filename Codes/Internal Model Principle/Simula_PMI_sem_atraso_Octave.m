clc
clear all
close all

%% Define caso (0 - nominal, 1 -caso incerto)
flag=0;

% Variável complexa s (função de transferência)
s=tf([1 0],1);

%% Sistema simulado
P=1/((s+1)*(0.5*s+1)*(0.25*s+1)*(0.125*s+1));
L=0.1; % atraso

%% Modelo nominal para fins de projeto
Pn=1/((1.15*s+1)*(0.8*s+1));
Ln=0; % atraso

P=Pn;

%% Tempo total de simulação 
tsim=50; % Tempo de simulação

%% Tempo de alteração da
tset=2;  % Instante de alteração da referência

%% Define tipo de perturbação (0 - degrau, 1-rampa, 2 senoide)

tipo_u=0;   % tipo de perturbação da entrada
tipo_y=0;   % tipo de perturbação da saída

tpu=25;     % Intante inicial da perturbação da entrada
du=inf;      % Intervalo de duração da perturbação da entrada
Apu=-0.2;   % Amplitude da perturbação da entrada

tpy=10;     % Intante inicial da perturbação da saída
dy=10;       % Duração da perturbação da saída
Apy=-0.2;    % Amplitude da perturbação da entrada

omega_u=2;  % Frequêcia da perturação da entrada (apenas para tipo 2)
omega_y=2;  % Frequêcia da perturação da saída   (apenas para o tipo 2)

%% Controlador 
kc=(1/0.3)^2;
C=kc*(1.15*s+1)*(0.8*s+1)/(s*(s+2/0.3));

%% Filtro de referência
F=(0.3*s+1)^2/(0.5*s+1)^2;

%% Define passo de integração e passo de armazenamonto 
h=0.001; % passo de integração
dT=0.1;  % passo de 

%%  Define passo automaticamete (0 - automatico, 1 - especificado)
flag_passo=1;

%%----------------- finaliza definição de parâmetros ------------------------

% Tolerância (ajuste de passo)
tol=1e-3;

if flag_passo==0;
 
  polo_real=max(real(pole(feedback(C*P))));
  dT=abs(polo_real)/20;
  h=dT/100;

end

if flag==0
  P=Pn;
  L=Ln;
end

% Representação em espaço de estados (acelera a simulação)
ssP=ss(P); Ap=ssP.a; Bp=ssP.b; Cp=ssP.c; Dp=ssP.d;
ssPn=ss(Pn); Apn=ssPn.a; Bpn=ssPn.b; Cpn=ssPn.c; Dpn=ssPn.d;
ssF=ss(F); Af=ssF.a; Bf=ssF.b; Cf=ssF.c; Df=ssF.d;
ssC=ss(C); Ac=ssC.a; Bc=ssC.b; Cc=ssC.c; Dc=ssC.d;


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
xf0=zeros(size(Af,1),1);
xc0=zeros(size(Ac,1),1);
xp0=zeros(size(Ap,1),1);

% Simula sistema de controle
 
for k=1:kfinal
  
   t(k)=(k-1)*dT;  % marca o tempo  

   % Estabelece as condições iniciais (apenas no primeiro passo)
   
   if k==1
    
    xp=xp0;   % condição inicial do processo
    xf=xf0;   % condição inicial do filtro
    xc=xc0;   % condição inicial do controlador     
    y=Cp*xp0;  % condição inicial da saída  (assumindo Dp=0)  
    
 
   end
 

 % Define mudança de referência
  if k<tset/dT+1
    r(k)=0;
  else
    r(k)=1; 
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
  
  % Implementa ajuste de passo
  
  if k==1
      
      mbackup=m;
      hbackup=h;
      
  elseif k>2
      
      if abs(Y(k-1)-Y(k-2))<tol          
          m=2;
          h=dT/2;          
      else
         m=mbackup;
         h=hbackup;
     end
     
 end
  
  
  
   for l=1:m
        
      
     %Calcula referência filtrada
    
     %[rf,xf]=sim_euler(Af,Bf,Cf,Df,xf,h,r(k));      
     
     xf=h*Af*xf+xf+h*Bf*r(k);
     rf=Cf*xf+Df*r(k);
     
     erro=rf-y; % erro com relação a referência filtrada
     
     %Calcula sinal de controle
     
    % [u,xc]=sim_euler(Ac,Bc,Cc,Dc,xc,h,erro);
     xc=h*Ac*xc+xc+h*Bc*erro;
     u=Cc*xc+Dc*erro;
    
     
       
   % [y_aux,xp]=sim_euler(Ap,Bp,Cp,Dp,xp,h,uatraso+qu);   
     xp=h*Ap*xp+xp+h*Bp*(u+qu);
     y_aux=Cp*xp+Dp*(u+qu);
     y=y_aux+qy;
     
     % Armazena valores calculados
     if l==1
       
       Xp(:,k)=xp;   
       Xf(:,k)=xf;  
       Xc(:,k)=xc;   
    
       U(k)=u;
       Y(k)=y;
       Qy(k)=qy;
       Qu(k)=qu;
    
    end
       
       
   end
   
    
   
end


nfont=15;   % define tamanho da fonte de texto
nlinha=3;   % define espessura da linha

figure(1)

subplot(2,1,1)
plot(t,r,'--r','linewidth',nlinha)
hold on
plot(t,Y,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(Y)-0.2  max(Y)+0.2])
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
plot(t,U,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(U)-0.2 max(U)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Controle');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)


figure(2)

subplot(2,1,1)
plot(t,Qy,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(Qy)-0.2  max(Qy)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Perturbação - Saída');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

subplot(2,1,2)
plot(t,Qu,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(Qu)-0.2 max(Qu)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Perturbação - Entrada');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

