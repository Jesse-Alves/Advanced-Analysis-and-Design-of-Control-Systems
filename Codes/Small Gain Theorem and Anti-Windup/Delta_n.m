function [G,Gn,Delta,w,s,Ln,Kn,tau_n] = Delta_n(L,K,tau)
    w=logspace(-2,4,10000);
    s=j*w;
    
    Ln = mean(L);
    Kn = mean(K);
    tau_n = mean(tau);
    
    Gn = Kn.*exp(Ln*s)./(tau_n*s+1);
    
    for i = 1:length(L)
        G(i,:) = K(i).*exp(-L(i)*s)./(tau(i)*s+1);
        Delta(i,:) = abs(1-G(i,:)./Gn);
    end
end    