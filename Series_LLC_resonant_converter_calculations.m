clc;

% Parameters
m = 8;
Q = 0.6;

Np = 1;          % Primary turns
Ns = 100;        % Secondary turns
Vo = 20e3;       % Output voltage in Volts
Po_max = 300;    % Max output power in Watts
Qmax = Q;        % Quality factor (already given)
fr = 50e3;       % Resonant frequency (Hz)

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


% Eq. 4: Calculate Rac_min
Rac_min = (8/pi^2) * (Np^2 / Ns^2) * (Vo^2 / Po_max);
fprintf('Rac_min = %.2f Ohms\n', Rac_min);

% Step 1: Solve intermediate values
A = (Qmax * Rac_min).^2;   % From Eq. 4 -> Lr / Cr = A
B = 1 / (2 * pi * fr).^2;  % From Eq. 5 -> Lr * Cr = B

% Step 2: Solve using symbolic substitution
% Lr / Cr = A => Lr = A * Cr
% (A * Cr) * Cr = B => A * Cr^2 = B => Cr = sqrt(B / A)
Cr = sqrt(B / A);
Lr = A * Cr;

% Display results
fprintf('Calculated Cr = %.3e F\n', Cr);
fprintf('Calculated Lr primary = %.3e H\n', Lr);

% Calculate Lm using m
Lm = Lr * (m - 1);
fprintf('Lm primary = %.2e H\n', Lm);