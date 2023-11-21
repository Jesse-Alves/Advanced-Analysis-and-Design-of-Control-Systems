%%SIMULADOR DE SISTEMA COM CONTROLE AMOSTRADO
%Aluno: Jess� Alves

%% Define caso (0 - nominal, 1 -caso incerto)
flag=0;

% Variavel s
s=tf([1 0],1);

%% Planta do Sistema
P=tf([-0.2 2],[1 2 1]);
L=0.0; % atraso na planta, se houver

%% Modelo nominal para fins de projeto
Pn=P;
Ln=0; % atraso

%% Instantes de Tempo
%Tempo total de simulacao
tsim=60; 

%Tempo da Mudanca de Referencia
tset=0;  

%% Tipo de perturbacao (0 - degrau, 1-rampa, 2 senoide)

tipo_u=0;   % tipo de perturbacao da entrada
tipo_y=0;   % tipo de perturbacao da saida

tpu=40;     % Intante inicial da perturbacao da entrada
du=inf;      % Intervalo de duracao da perturbacao da entrada
Apu=0.2;   % Amplitude da perturbacao da entrada

tpy=20;     % Intante inicial da perturbacao da saida
dy=inf;       % Duracao da perturbacao da saida
Apy=0.2;    % Amplitude da perturbacao da entrada

omega_u=2;  % Frequencia da perturacao da entrada (apenas para tipo 2)
omega_y=2;  % Frequecia da perturacao da saida   (apenas para o tipo 2)

%% Controlador 

% Periodo de amostragem
Ts=ta; 
% Modelo discretizado
Pd=c2d(P,Ts,'zoh');                       

tipo_cont = 0;
if tipo_cont == 0
    % =========== Controle por Alocacao de Polos ===========
    C = Cz;
    F = Fz;
    Rp = denCz;
    Sp = numCz;
    Tp = numFz;
else
    % ==================== Controle IMC ====================
    C = C_IMC; 
    [numC_IMC denC_IMC] = tfdata(C,'v');
    %F=1;    
    Rp = denC_IMC;
    Sp = numC_IMC;
    Tp = Sp;
end

Hmf = minreal(feedback(C*Pd,1));

% Torna o polinomio monico
Rp=Rp/Rp(1);
Sp=Sp/Rp(1);
Tp=Tp/Rp(1);

% Polinomio modificado para considerar controles passados
Rm=Rp(2:size(Rp,2));
nr=size(Rp,2)-1;
ns=size(Sp,2);
nt=size(Tp,2);

%% Define passo de integracao e passo de armazenamento 
dT=Ts;     % passo de amostragem 
h=Ts/100;  % passo de integracao

%%----------------- finaliza definicao de parametros ------------------------
if flag==0
  P=Pn;
  L=Ln;
end

% Representação em espaço de estados (acelera a simulação)
ssP=ss(P); Ap=ssP.a; Bp=ssP.b; Cp=ssP.c; Dp=ssP.d;
ssPn=ss(Pn); Apn=ssPn.a; Bpn=ssPn.b; Cpn=ssPn.c; Dpn=ssPn.d;

% Arredonda o passo de integração tal que dT seja múltiplo de h
m=ceil(dT/h);
h=dT/m;
% Atraso de tempo discretro
Ld=round(L/h);
U_atraso=zeros(1,Ld+1);

U_passado=zeros(1,nr);
Yr_passado=zeros(1,nt);
Y_passado=zeros(1,ns);

%display('Atraso multiplo do passo - atraso considerado:')
%display(Ld*h);

% Calcula numero de passos de simulacao
kfinal=ceil(tsim/dT);

% Estabelece condicoes iniciais dos elementos dinamicos
xp0=zeros(size(Ap,1),1);

% Simula sistema de controle 
for k=1:kfinal
  
   t(k)=(k-1)*dT;  % marca o tempo  

   % Estabelece as condições iniciais (apenas no primeiro passo)
   
   if k==1    
    xp=xp0;   % condição inicial do processo
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

   Yr_passado=[r(k) Yr_passado(1:nt-1)];
   
   % R(q^-1)u(k)=T(q^-1)y_r(k)-S(q^-1)y(k)
    u=-Rm*U_passado'+Tp*Yr_passado'-Sp*Y_passado';
  
   for l=1:m

     %Calcula sinal de controle    

     %Calcula saida do processo
     U_atraso=[u U_atraso(1:Ld)];
     u_atraso=U_atraso(Ld+1);
     
   % [y_aux,xp]=sim_euler(Ap,Bp,Cp,Dp,xp,h,uatraso+qu);   
     xp=h*Ap*xp+xp+h*Bp*(u_atraso+qu);
     y_aux=Cp*xp+Dp*(u_atraso+qu);
     y=y_aux+qy;

     % Armazena valores calculados
     if l==1       
       Xp(:,k)=xp;            
       U(k)=u_atraso;
       Y(k)=y;
       Qy(k)=qy;
       Qu(k)=qu;    
    end         
   end
     Y_passado=[y Y_passado(1:ns-1)];
     U_passado=[u U_passado(1:nr-1)];
end

nfont=15;   % define tamanho da fonte de texto
nlinha=3;   % define espessura da linha

%% Plotagem dos graficos 
% ===> Sinal de Saida e Referencia
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

% ===> Sinal de Controle
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

% ===> Sinais de Pertubacao
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


