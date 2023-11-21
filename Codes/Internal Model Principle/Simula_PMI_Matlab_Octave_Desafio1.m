clc
clear all
close all

%% Define caso (0 - nominal, 1 -caso incerto)
flag=1;

% Variável complexa s (função de transferência)
s=tf([1 0],1);

%% Sistema simulado
P=2/(s+1.5); % <==== Muda Aqui
L=0; % atraso

%% Modelo nominal para fins de projeto
Pn=1/((1.15*s+1)*(0.8*s+1));
Ln=0; % atraso



%% Tempo total de simulação 
tsim=60; % Tempo de simulação

%% Tempo de alteração da
tset=2;  % Instante de alteração da referência

%% Define tipo de perturbação (0 - degrau, 1-rampa, 2 senoide)

tipo_u=0;   % tipo de perturba��o da entrada
tipo_y=0;   % tipo de perturba��o da sa�da

tpu=40;     % Instante inicial da perturba��o da entrada
du=inf;      % Intervalo de duração da perturba��o da entrada
Apu=-0.2;   % Amplitude da perturba��o da entrada

tpy=20;     % Intante inicial da perturba��o da sa�da
dy=inf;       % Duração da perturba��o da sa�da
Apy=-0.2;    % Amplitude da perturba��o da sa�da

omega_u=2;  % Frequêcia da perturba��o da entrada (apenas para tipo 2)
omega_y=2;  % Frequêcia da perturba��o da saída   (apenas para o tipo 2)

%% Controlador 
kc=1.25;
z = 1.35;
C=(kc*(s+z))/(s);

%% Filtro de referência
F=1;

infor = stepinfo(series(feedback(C*P,1),F))
polos = pole(series(feedback(C*P,1),F))



%% Define passo de integração e passo de armazenamonto 
h=0.001; % passo de integração
dT=0.1;  % passo de 

%%  Define passo automaticamete (0 - automatico, 1 - especificado)
flag_passo=0;

%%----------------- finaliza definição de parâmetros ------------------------

% Tolerância (ajuste de passo)
tol=1e-3;

if flag_passo==0;
 
  polo_real=max(real(pole(feedback(C*P,1))));
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
display('Atraso m�ltiplo do passo - atraso considerado:')
display(Ld*h)

% Calcula n�mero de passos de simula��o
kfinal=ceil(tsim/dT);

% Estabelece condi��es iniciais dos elementos din�micos
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
    
     
     %Calcula saída do processo
     U_atraso=[u U_atraso(1:Ld)];
     u_atraso=U_atraso(Ld+1);
     
   % [y_aux,xp]=sim_euler(Ap,Bp,Cp,Dp,xp,h,uatraso+qu);   
     xp=h*Ap*xp+xp+h*Bp*(u_atraso+qu);
     y_aux=Cp*xp+Dp*(u_atraso+qu);
     y=y_aux+qy;
     
     % Armazena valores calculados
     if l==1
       
       Xp(:,k)=xp;   
       Xf(:,k)=xf;  
       Xc(:,k)=xc;   
    
       U(k)=u_atraso;
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
plot(t,r,'--r','linewidth',nlinha-0.5)
hold on
plot(t,Y,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(Y)-0.2  max(Y)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Sa�da');
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
hy=ylabel('Perturba��o - Sa�da');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

subplot(2,1,2)
plot(t,Qu,'b','linewidth',nlinha)
grid on
axis([0 t(end) min(Qu)-0.2 max(Qu)+0.2])
set(gca,'fontsize',nfont)

% Define as entradas de texto e ajusta o tamanho da fonte
hx=xlabel('Tempo (s)');
hy=ylabel('Perturba��o - Entrada');
set(hx,'fontsize',nfont)
set(hy,'fontsize',nfont)

