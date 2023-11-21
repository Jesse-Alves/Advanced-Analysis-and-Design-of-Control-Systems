function [xkplus1, y] = get_ss_output(xk, ss_matrix = ss(0,0,0,0), uk)
xkplus1 = (ss_matrix.a*xk) + (ss_matrix.b*uk);
y = ss_matrix.c*xkplus1 + (ss_matrix.d*uk);
end


% Considering Y(z) = H(z)X(z) -> H(z)=Y(z)/X(z)=b(z)/a(z)
% y[k] = ([b0 b1 ... bm]*(x[k] x[k-1] ... x[k-m])
%       - [a1 b2 ... an]*(y[k-1] y[k-2] ... y[k-n]))/a0
function ykplus1 = compute_difference_equation_step(y, x, a, b)
    ykplus1 = (b*x' - a(2:end)*y')/a(1);
endfunction


% Discrete simulation based in SANTOS, T. L. M. code example 
function [U, Y, E, R] = simulate_discrete_sys(sim_time, dt, integration_step_ratio,
                                          plant_tf, digital_controller,
                                          reference_filter, reference,
                                          input_disturbance = 0,
                                          output_disturbance = 0)

has_filter = false;
if length(pole(reference_filter)) > 0
    has_filter = true;
    [reference_filter_num, reference_filter_den, ts] = filtdata(reference_filter);
    [reference_filter_num, reference_filter_den] = deal(reference_filter_num{1},
                                                        reference_filter_den{1})
    reference_filtered = zeros(1, (length(reference_filter_den) - 1));   % Initial conditions
    ref_filt_num_delays = length(reference_filter_num)-1;
    ref_filt_den_delays = length(reference_filter_den)-1;
    reference = [zeros(1, length(reference_filter_num) - 1) reference']; % Add r[-1], r[-2], r[-n-1]
end


plant_ss = ss(plant_tf);
plant_ss.a = eye(size(plant_ss.a)) + dt*plant_ss.a;
plant_ss.b = dt*plant_ss.b;
plant_states = [zeros(size(plant_ss.a,1),1)];
Y(1) = plant_ss.c*plant_states;


[digital_controller_num, digital_controller_den, sampling_period] = filtdata(
    digital_controller);
[digital_controller_num, digital_controller_den] = deal(
    digital_controller_num{1}, digital_controller_den{1});
dc_num_delays = length(digital_controller_num)-1;
dc_den_delays = length(digital_controller_den)-1;
digital_controller_out = zeros(1, (length(digital_controller_den) - 1));   % Initial conditions
err = zeros(1, length(digital_controller_num)-1);   % Initial conditions


for k=1:length(sim_time)
    if (mod(k, integration_step_ratio+1) == 0) || (k == 1)
        digital_controller_out_idx = length(digital_controller_out) + 1;
        err_idx = length(err) + 1;
        if !has_filter
            err(err_idx) = reference(k) - Y(k);
        else
            rf_k = length(reference_filtered)+1;
            reference_filtered(rf_k) = compute_difference_equation_step(
                                                reference_filtered(end:-1:end-(ref_filt_den_delays-1)),
                                                reference(k+ref_filt_num_delays:-1:k),
                                                reference_filter_den, 
                                                reference_filter_num);
            err(err_idx) = reference_filtered(rf_k) - Y(k);
        end
        digital_controller_out(digital_controller_out_idx) = compute_difference_equation_step(
            digital_controller_out(end:-1:end-(dc_den_delays-1)),
            err(end:-1:end-dc_num_delays), digital_controller_den,
            digital_controller_num);
    end


    U(k) = digital_controller_out(digital_controller_out_idx) + input_disturbance(k);
    [plant_states, y_aux] = get_ss_output(plant_states, plant_ss, U(k));
    Y(k+1) = y_aux + output_disturbance(k);
    E(k) = err(err_idx);
    if has_filter
        R(k) = reference_filtered(rf_k);
    else
        R = reference;
    end
    end
Y = Y(2:end); % Remove initial condition y(0)
end