clc;

% Parameters
m = 8; % previous value = 8
Q = 0.105; % previous val = 0.6

Np = 1;          % Primary turns
Ns = 75;        % Secondary turns
Vo = 15*1e3;       % Output voltage in Volts
Po_max = 300;    % Max output power in Watts
Qmax = Q;        % Quality factor (already given)
fr = 20e3;       % Resonant frequency (Hz)
R0 = 30e3;
Cr = 6.8 * 1e-6;

Fx = linspace(0.1, 10.0, 10000); % Sweep Fx from 0 to 50

% Eq. 4: Calculate Rac_min
Rac_min = (8/pi^2) * (Np^2 / Ns^2) * (Vo^2 / Po_max);
%Rac_min = (8/pi^2) * (Np^2 / Ns^2) * R0;
fprintf('Rac_min = %.2f Ohms\n', Rac_min);


% Step 1: Solve intermediate values
%A = (Qmax * Rac_min).^2;   % From Eq. 4 -> Lr / Cr = A
%B = 1 / (2 * pi * fr).^2;  % From Eq. 5 -> Lr * Cr = B

% Step 2: Solve using symbolic substitution
% Lr / Cr = A => Lr = A * Cr
% (A * Cr) * Cr = B => A * Cr^2 = B => Cr = sqrt(B / A)
%Cr = sqrt(B / A);
%Lr = A * Cr;
Lr = ((1 / (2 * pi * fr)).^2) / Cr;

%Qmax = sqrt(Lr / Cr) / Rac_min;
%Q = Qmax;


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


% Display results
fprintf('Calculated Cr = %.3e F\n', Cr);
fprintf('Calculated Lr primary = %.3e H\n', Lr);

% Calculate Lm using m
Lm = Lr * (m - 1);
fprintf('Lm primary = %.2e H\n', Lm);

% Calculate turns ratio
a = Np / Ns;

% Calculate secondary side inductances
Lm_secondary = Lm / (a^2);
Llk_secondary = Lr / (a^2);

% Display results
fprintf('Secondary Magnetizing Inductance: %.6e H\n', Lm_secondary);
fprintf('Secondary Leakage Inductance:     %.6e H\n', Llk_secondary);


fprintf('Qmax=%.4f\n\n',Q);

fprintf('.param Cr=%.4fu\n',Cr*1e6);
fprintf('.param Lr=%.4fu\n',Lr*1e6);
fprintf('.param Lm=%.4fu\n',Lm*1e6);
fprintf('.param LSm=%.4fm\n',Lm_secondary*1e3);
fprintf('.param LSr=%.4fm\n',Llk_secondary*1e3);
