% MERTONWEAK  Test weak convergence of EM and Milstein for PD
% Based closely on Higham  emweak.m

clear all; close all; rng(100);

% Merton Parameters
mu = 0.05; sigma = 0.2; V0 = 100; T = 1; K = 80; 
M = 9000000; % Increased M
p_max = 5;

% Analytical PD (Ground Truth)
d2 = (log(V0/K) + (mu - 0.5*sigma^2)*T) / (sigma*sqrt(T));
true_PD = 0.5 * (1 + erf(-d2 / sqrt(2)));

Dt_vals = 2.^- (1:p_max);
err_em = zeros(p_max,1);
err_mil = zeros(p_max,1);

for p = 1:p_max
    Dt = Dt_vals(p);
    N = round(T/Dt);
    V_em = V0*ones(M,1);
    V_mil = V0*ones(M,1);
    
    for n = 1:N
        dW = sqrt(Dt)*randn(M,1);
        % Vectorized updates
        V_em = V_em + mu*V_em*Dt + sigma*V_em.*dW;
        V_mil = V_mil + mu*V_mil*Dt + sigma*V_mil.*dW + ...
                0.5*sigma^2*V_mil.*(dW.^2 - Dt);
    end
    
    % Weak Error: Absolute difference in PD
    err_em(p) = abs(true_PD - mean(V_em < K));
    err_mil(p) = abs(true_PD - mean(V_mil < K));
end

% Plotting as per Higham's Fig 5.1
figure;
loglog(Dt_vals, err_em, 'b-o', Dt_vals, err_mil, 'r-x'), hold on
loglog(Dt_vals, Dt_vals, 'k--') % Reference Slope 1.0
set(gca, 'FontSize', 28);
xlabel('\Delta t', 'FontSize', 28); 
ylabel('Weak Error', 'FontSize', 28);
legend('EM', 'Milstein', 'Order 1.0', 'Location', 'Best', 'FontSize', 26);
grid on