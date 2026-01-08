% MERTONSTRONG  Test strong convergence of EM and Milstein
% Based on Higham scripts emstrong.m and milstrong.m

clear all; rng(100);
mu = 0.05; sigma = 0.2; V0 = 100; T = 1; % Merton Parameters
M = 1000;                                % number of paths
p_max = 6;                               % max power of 2 for dt

Dt_vals = 2.^-(1:p_max);                 % range of timesteps
err_em = zeros(p_max,1);
err_mil = zeros(p_max,1);

for p = 1:p_max
    Dt = Dt_vals(p);
    N = T/Dt;
    V_em = V0*ones(M,1);
    V_mil = V0*ones(M,1);
    W_T = zeros(M,1);
    
    for n = 1:N
        dW = sqrt(Dt)*randn(M,1);
        W_T = W_T + dW;
        
        % Euler-Maruyama
        V_em = V_em + mu*V_em*Dt + sigma*V_em.*dW;
        
        % Milstein (using g(V) = sigma*V, so g'g = sigma^2*V)
        V_mil = V_mil + mu*V_mil*Dt + sigma*V_mil.*dW ...
                + 0.5*sigma^2*V_mil.*(dW.^2 - Dt);
    end
    
    % Exact solution at T (Ground Truth)
    V_exact = V0*exp((mu - 0.5*sigma^2)*T + sigma*W_T);
    
    err_em(p) = mean(abs(V_exact - V_em));
    err_mil(p) = mean(abs(V_exact - V_mil));
end

% Plotting as per Higham's Fig 4.2
figure;
loglog(Dt_vals, err_em, 'b-o', Dt_vals, err_mil, 'r-x'), hold on
loglog(Dt_vals, Dt_vals.^0.5, 'k--'), loglog(Dt_vals, Dt_vals, 'k--')
set(gca, 'FontSize', 28);
xlabel('\Delta t', 'FontSize' , 28), 
ylabel('Strong Error', 'FontSize',28)
legend('EM','Milstein','Order 0.5','Order 1.0','Location','Best','FontSize',28)
grid on