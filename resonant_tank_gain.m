clc;

% Parameters
m = 8;
Q = 0.6;
Fx = linspace(0, 10, 1000); % Sweep Fx from 0 to 10

% Equation from image
numerator = Fx.^2 .* (m - 1);
denominator = sqrt((m .* Fx.^2 - 1).^2 + Fx.^2 .* (Fx.^2 - 1).^2 .* (m - 1).^2 .* Q^2);
K = numerator ./ denominator;

% Find maximum value of K and corresponding Fx
[K_max, idx_max] = max(K);
Fx_max = Fx(idx_max);

% Plotting
figure;
plot(Fx, K, 'LineWidth', 2);
hold on;
plot(Fx_max, K_max, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Mark max point
grid on;
xlabel('F_x');
ylabel('K(Q, m, F_x)');
title('Plot of K(Q, m, F_x) vs F_x');
legend('K(F_x)', sprintf('Max: K=%.3f at F_x=%.3f', K_max, Fx_max), 'Location', 'best');

% Print results
fprintf('Maximum K = %.5f occurs at F_x = %.5f\n', K_max, Fx_max);