% Controlador
%K_lqr = [0,0]; para u=0
%ueq = 0; para u=0
%xeq = [0;0]; para u=0

%xeq = [1;0]; %para modelo 1
%ueq = 1; %para modelo 1
%K_lqr = [2.3166,3.8253]; %para modelo 1

%xeq = [4;0]; %para modelo 2
%ueq = 4; %para modelo 2
%K_lqr = [2.3166,1.4062]; %para modelo 2

% Euler

h = 0.01; % Passo de integração
Tmax = 20; % Tempo de simulação
t = 0:h:Tmax;
kfinal = size(t,2);
x = zeros(2,kfinal); 
u = zeros(1,kfinal);


x0 = [-5;-5]; % Condicao inicial


    % Trajetória do oscilador

x(:,1) = x0;
for k = 1:kfinal-1
    u(k) = -K_lqr*(x(:,k)-xeq)+ueq;
    dx = h*[x(2,k); -x(1,k)+0.3*(1-x(1,k)^2)*x(2,k)+u(k)];
    x(:,k+1) = x(:,k) + dx;
end

% Plot da resposta 
figure(1)
plot(t,x(1,:),'r','linewidth',1);
hold on
plot(t,x(2,:),'b','linewidth',1);
plot(t,xeq(1)*ones(1,kfinal),'m--','linewidth',1)
xlabel('Tempo (s)')
ylabel('Estados do sistema')
legend({'x_1(t)','x_2(t)'}, 'FontSize',8)