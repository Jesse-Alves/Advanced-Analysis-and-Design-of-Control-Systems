



xk0 = zeros(size(xeq_continue));
a_euler = eye(size(gs_ss.A))+(integration_step_size*gs_ss.A);
b_euler = integration_step_size*gs_ss.B;
xt(:,1) = xk0;
xd(:,1) = xk0;
ut(:,1) = u0;

k=2;

integration_step_ratio = ta/


for i=1:sim_time_length
    xt(:,i+1) = a_euler*xt(:,i) + b_euler*ut(:,i);
    yt(:,i) = gs_ss.c*xt(:,i);

    % Momento de amostragem do sistema incluido o tempo zero
    if (mod(i, integration_step_ratio) == 1) || i == 1
        % Apenar para verificar que amostragem é realizada nos momentos
        % corretos. Para verificar, basta remover o ponto e vírgula
        
        t1 = (i-1)*integration_step_size;
        t2 = (k-1)*sampling_period;
        
        xd(:,k) = xt(:,i);
        yd(:,k) = yt(:,i);
        ud(:,k) = -kd*(xd(:,k) - xeq_discrete.*unit_step(:,i)) + ueq(:,i);

        k = k+1;
    endif

    if i != sim_time_length
        ut(:,i+1) = ud(:,k-1);
    endif
endfor
