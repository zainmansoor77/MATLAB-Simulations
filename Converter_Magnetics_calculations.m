clc;

B = 4;      % flux density (B) in Tesla
B_sat = 1;  % saturation flux density (B_sat) in Tesla
Ac = 1;     % cross-sectional area (A_c) in m²
U0 = 4 * pi * 1e-7; % Permeability of free space [H/m]
Ur = 1;     % relative permeability 
i = 1;      % winding current (i) in Amperes
lm = 1;     % core length (l_m) in meters
n = 1;      % number of turns

% Calculate magnetic flux (Phi)
Phi = B * Ac;

% Compute magnetic field strength (H)
H = (n * i) / lm; % [A/m]

% Linear model (no saturation)
U = Ur * U0; % Core permeability [H/m]
B_linear = U0 * Ur * H;

% Piecewise model (with saturation)
if abs(B_linear) <= B_sat
    B_piecewise = B_linear;
else
    B_piecewise = sign(B_linear) * B_sat; % Clamp to B_sat
end

%% Derived Parameters
I_sat = (B_sat * lm) / (U * n); % Saturation current [A] 
L = (U * n^2 * Ac) / lm; % Inductance [H]


%%  Display flux result
fprintf('\nMagnetic Flux (Φ): %.4f Weber (Wb)\n', Phi);
%% Display Results for Input Current
fprintf('\n--- Results for Input Current ---\n');
fprintf('Magnetic Field Strength (H): %.2f A/m\n', H);
fprintf('Linear Model Flux Density (B_linear): %.9f T (%.3f uT) \n', B_linear, B_linear * 1e6);
fprintf('Piecewise Model Flux Density (B_piecewise): %.9f (%.3f uT) T\n', B_piecewise, B_piecewise * 1e6);

%% Plot B-H Characteristics
% Generate H values for plotting
H_max = (B_sat / (U0 * Ur)) * 1.5; % Extend beyond saturation
H_range = linspace(0, H_max, 1000); % Range of H [A/m]

% Linear model (B = μ₀μ_r*H)
B_linear_range = U0 * Ur * H_range;

% Piecewise model (saturation)
B_piecewise_range = U0 * Ur * H_range;
B_piecewise_range(abs(B_piecewise_range) > B_sat) = ...
    sign(B_piecewise_range(abs(B_piecewise_range) > B_sat)) * B_sat;

% Plot
figure;
plot(H_range, B_linear_range, 'b--', 'LineWidth', 1.5); hold on;
plot(H_range, B_piecewise_range, 'r', 'LineWidth', 1.5);
xline(H, 'k:', 'LineWidth', 1.5); % Mark input H value
yline(B_sat, 'k--', 'LineWidth', 1.5); % Mark B_sat
xlabel('Magnetic Field Strength (H) [A/m]');
ylabel('Flux Density (B) [T]');
title('B-H Characteristic of Magnetic Core');
legend('Linear Model', 'Piecewise Model (Saturation)', 'Input H', 'B_{sat}', 'Location', 'southeast');
grid on;

%% Display Saturation current and Inductance Parameters
fprintf('\n--- Key Parameters ---\n');
fprintf('Inductance (L): %.9f H (%.3f uH)\n', L, L * 1e6);
fprintf('Saturation Current (I_sat): %.4f A\n', I_sat);