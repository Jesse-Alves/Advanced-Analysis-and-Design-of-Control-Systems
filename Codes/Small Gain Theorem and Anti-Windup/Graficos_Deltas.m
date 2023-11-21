L = [0.9 0.7 0.6 0.4];
K = [1.3 0.9 1.2 0.8];
tau = [1.2 1.1 0.8 0.9];
tau_c = 0.5;

[G,Gn,Delta,w,s,Ln,Kn,tau_n] = Delta_n(L,K,tau);
for i = 1:length(L)
    titulo = 'Modulo da incerteza multiplicativa ' + i;
    figure;
    semilogx(w,Delta(i,:));
    title(titulo);
    xlabel('Frequência w','Fontsize',14);
    ylabel('Incerteza Multiplicativa','Fontsize',14);
    grid on
end

Kc = tau_n/(Kn*(tau_c + Ln));
Ti = min(tau_n,4*(tau_c + Ln));
Cs = Kc*(s*Ti + 1)./s*Ti;

Comp_s = (Gn.*Cs)./(1+Gn.*Cs);
Delta_b = max(Delta,[],1);

figure;
semilogx(w,abs(Comp_s).*Delta_b);
title('Função complementar de sensibilidade multiplicada pelas incertezas');
xlabel('Frequência w','Fontsize',14);
ylabel('|C(jw)|Delta-barra(w)','Fontsize',14);
grid on;


