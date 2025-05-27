clc;

B = 4;      % flux density (B) in Tesla
B_sat = 1;  % saturation flux density (B_sat) in Tesla
Ac = 1;     % cross-sectional core area (Ac) in m²
Ag = 0.1;   % cross-sectional core air gap area (Ag) in m²
U0 = 4 * pi * 1e-7; % Permeability of free space [H/m]
Ur = 1;     % relative permeability 
i = 1;      % winding current (i) in Amperes
I_rms = 1;  % RMS current
lm = 1;     % core length (lm) in meters
lg = 0.1;   % air gap length (lg) in meters
n = 3;      % number of turns
n1 = 1;     % primary turns (n1)
n2 = 1;     % secondary turns (n2)
i1 = 1;     % primary current (i1) in A
i2 = 1;     % secondary current (i2) in A
L_M = 1;    % Primary Magnetizing inductance
Ll1 = 1;    % Primary Leakage Inductance
Ll2 = 1;    % Secondary Leakage Inductance
L12_choice = 1; % 1 for  L12 = (n1 * n2) / R; 2 for L12 = (n2 / n1) * L_M;
K_fe = 1;   % material constant
beta = 1;   % material constant
delta_B = 1;% peak ac flux density (ΔB) in Tesla
MLT = 0.15;             % Mean Length per Turn [m]
AWG = 24;               % American Wire Gauge [] % diameter of the conductor is calculated with the AWG
h = 0.5;                % height (or thickness) of the conductor layer in the direction perpendicular to the surface where the current flows
temp = 100;             % Temperature [°C]
rho_20C = 1.724e-8;     % Copper resistivity at 20°C [Ω-m]
f = 100e3;              % Frequency [Hz]


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

%%  Reluctances Calculations
% 1. Compute reluctances
Rc = lm / (Ur * U0 * Ac); % Core reluctance (Eq. 10.25)
Rg = lg / (U0 * Ag); % Air gap reluctance (Eq. 10.25)
R_total = Rc + Rg; % Total reluctance (series) with air gap
R_total_no_gap = Rc; % Total reluctance without air gap

% 2. Compute total MMF (source)
MMF_total = n * i; % Total MMF (nI) [A-turns]

% 3. Compute magnetic flux (Phi)
Phi = MMF_total / R_total; 

% 4. Compute MMF drops across core and air gap
MMF_core = Phi * Rc; % MMF across core
MMF_airgap = Phi * Rg; % MMF across air gap

% 5. Validate MMF_total = MMF_core + MMF_airgap
MMF_error = abs(MMF_total - (MMF_core + MMF_airgap));

% 2. Inductance (Eq. 103)
L_with_gap = n^2 / R_total; % Inductance with air gap
L_no_gap = n^2 / R_total_no_gap; % Inductance without air gap

% 3. Saturation current (Eq. 105)
I_sat_with_gap = (B_sat * Ac / n) * (Rc + Rg); % With air gap
I_sat_no_gap = (B_sat * Ac / n) * Rc; % Without air gap

% 4. Flux vs. MMF (ni) characteristics
ni = linspace(0, 2*I_sat_with_gap*n, 1000); % MMF range [A-turns]
Phi_linear_with_gap = ni / R_total; % Flux (linear region)
Phi_linear_no_gap = ni / R_total_no_gap; % Flux (linear region)

% Clamp flux at saturation (Phi_sat = B_sat * A_c)
Phi_sat = B_sat * Ac; % Eq. 104
Phi_with_gap = min(Phi_linear_with_gap, Phi_sat);
Phi_no_gap = min(Phi_linear_no_gap, Phi_sat);

%% Calculations for Real Transformer
    % 1. Compute core reluctance (Eq. 10.34)
    R = lm / (Ur * U0 * Ac); 

    % 2. Compute total MMF (Eq. 10.35)
    MMF = n1 * i1 + n2 * i2;

    % 3. Compute flux (Eq. 10.36)
    Phi = MMF / R;

% 2. Compute magnetizing inductance (Eq. 10.44)
L_M = n1^2 / R;

% 3. Compute magnetizing current (Eq. 10.44)
i_M = i1 + (n2/n1) * i2;

if L12_choice == 1
    L12 = (n1 * n2) / R; % Eq. 10.49 (using reluctance)
elseif L12_choice == 2
     L12 = (n2 / n1) * L_M; % Eq. 10.49 (using L_M)
end

% Self-inductances (Eq. 10.50)
L11 = Ll1 + (n1 / n2) * L12;
L22 = Ll2 + (n2 / n1) * L12;

% Effective turns ratio (Eq. 10.51)
n_e = sqrt(L22 / L11);

% Coupling coefficient (Eq. 10.52)
k = L12 / sqrt(L11 * L22);

% Calculate core loss
    P_fe = K_fe * (delta_B)^beta * Ac * lm; % Eq. 10.57

%% Wire Dimensions
% AWG to diameter conversion (from table)
dia = (0.127*92.^((36-AWG)/39))/1000; % Diameter in meters
A_w = pi*(dia/2)^2;     % Wire cross-sectional area [m²]

% Temperature Adjustment for Resistivity
rho = rho_20C*(1 + 0.00393*(temp - 20)); % Temperature coefficient for copper

% DC Resistance Calculation
R_dc = rho * (n * MLT) / A_w;
P_cu = I_rms^2 * R_dc;

% Skin Effect Analysis
skin_depth = sqrt(rho/(pi*f*U0*Ur)); % Standard skin depth formula

% Calculate AC resistance (Equation 10.63)
R_ac = (h / skin_depth) * R_dc;

% Calculate copper loss in layer 1 (Equation 10.64)
P_1 = i^2 * R_ac;

% Pre-calculate Common Factor
common_factor = I_rms^2 * (h/skin_depth) * R_dc;

% Calculate Loss for Each Layer
layer_copper_losses = zeros(1, n);
for m = 1:n
    layer_copper_losses(m) = common_factor * ((m-1)^2 + m^2);
end

% Calculate Total Copper Loss
P_cu1 = sum(layer_copper_losses);

% Equation Implementation
term1 = I_rms^2;                    % Current squared term
term2 = (h / skin_depth) * R_dc;            % Geometry/resistance factor
term3 = n/3;                    % Layer scaling factor
term4 = (2*n^2 + 1);            % Proximity effect polynomial

P_cu2 = term1 * term2 * term3 * term4;

P_dc = (I_rms^2) * n * R_dc;

Fr = (1/3) * (h/skin_depth) * (2*n^2 + 1);


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

%% Reluctance Results Display
fprintf('\n--- Reluctance Results Display ---\n');
fprintf('Core Reluctance (R_core): %.4e A-turns/Wb\n', Rc);
fprintf('Air Gap Reluctance (R_airgap): %.4e A-turns/Wb\n', Rg);
fprintf('Total Reluctance (R_total): %.4e A-turns/Wb\n', R_total);
fprintf('Total MMF (nI): %.2f A-turns\n', MMF_total);
fprintf('Magnetic Flux (Φ): %.4e Wb\n', Phi);
fprintf('MMF Core: %.2f A-turns\n', MMF_core);
fprintf('MMF Air Gap: %.2f A-turns\n', MMF_airgap);
fprintf('MMF Validation Error: %.2e A-turns\n', MMF_error);

%% Plot Reluctance Contributions (Optional)
figure;
bar([Rc, Rg]);
xticklabels({'Core', 'Air Gap'});
ylabel('Reluctance (A-turns/Wb)');
title('Reluctance Distribution');
grid on;

%% Display Results
fprintf('\n--- Results ---\n');
fprintf('Core Reluctance (R_c): %.4e A-turns/Wb\n', Rc);
fprintf('Air Gap Reluctance (R_g): %.4e A-turns/Wb\n', Rg);
fprintf('Inductance with Air Gap (L): %.4f H\n', L_with_gap);
fprintf('Inductance without Air Gap (L): %.4f H\n', L_no_gap);
fprintf('Saturation Current with Air Gap (I_sat): %.4f A\n', I_sat_with_gap);
fprintf('Saturation Current without Air Gap (I_sat): %.4f A\n', I_sat_no_gap);

%% Plot Results
figure;

% Subplot 1: Inductance and Saturation Current Comparison
subplot(2,1,1);
bar([L_no_gap, L_with_gap; I_sat_no_gap, I_sat_with_gap]);
xticks(1:2);
xticklabels({'Inductance (H)', 'Saturation Current (A)'});
legend('No Air Gap', 'With Air Gap', 'Location', 'northwest');
title('Effect of Air Gap on Inductance and Saturation Current');
grid on;

% Subplot 2: Flux vs. MMF (ni)
subplot(2,1,2);
plot(ni, Phi_no_gap, 'r--', 'LineWidth', 1.5); hold on;
plot(ni, Phi_with_gap, 'b', 'LineWidth', 1.5);
yline(Phi_sat, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Φ_{sat}');
xlabel('MMF (ni) [A-turns]');
ylabel('Flux (Φ) [Wb]');
title('Φ vs. MMF (ni)');
legend('No Air Gap', 'With Air Gap', 'Φ_{sat}', 'Location', 'southeast');
grid on;

 %% Display Transformer Results
    fprintf('\n--- Real Transformer Results ---\n');
    fprintf('Core Reluctance (R): %.4e A-turns/Wb\n', R);
    fprintf('Total MMF (n1i1 + n2i2): %.2f A-turns\n', MMF);
    fprintf('Magnetic Flux (Φ): %.4e Wb\n', Phi);
    fprintf('Magnetizing Inductance (L_M): %.4f H\n', L_M);
    fprintf('Magnetizing current (i_M): %.4f A\n', i_M);

%% Display Results
fprintf('\n--- Results ---\n');
fprintf('Mutual Inductance (L12): %.9f H (%.3f uH)\n', L12, L12 * 1e6);
fprintf('Primary Self-Inductance (L11): %.9f H (%.3f uH)\n', L11, L11 * 1e6);
fprintf('Secondary Self-Inductance (L22): %.9f H (%.3f uH)\n', L22, L22 * 1e6);
fprintf('Effective Turns Ratio (n_e): %.4f\n', n_e);
fprintf('Coupling Coefficient (k): %.4f\n', k);
if k >= 0.99
    fprintf('Note: High coupling (k ≈ 1), effective turns ratio ≈ physical ratio (n2/n1 = %.2f)\n', n2/n1);
end

 % Display result
    fprintf('\n--- Empirical Core Loss ---\n');
    fprintf('Total Core Loss (P_fe): %.4f W\n', P_fe);

%% Display Transforner Results
fprintf('\nTransformer Winding Loss Analysis:\n');
fprintf('---------------------------------\n');
fprintf('DC Resistance: %.3f mΩ\n', R_dc*1000);
fprintf('DC Copper Loss: %.2f W\n', P_cu);
fprintf('Skin Depth: %.3f mm\n', skin_depth*1000);
fprintf('AC Resistance (R_ac): %.3f mΩ\n', R_ac*1000);
fprintf('Copper Loss in Layer 1 (P_1): %.4f W\n', P_1);

% Display Results
fprintf('Multi-Layer Winding Copper Loss Analysis\n');
fprintf('----------------------------------------\n');
fprintf('Input Parameters:\n');
fprintf('RMS Current (I): %.1f A\n', I_rms);
fprintf('DC Resistance/layer (R_de): %.4f Ω\n', R_dc);
fprintf('Number of Layers (M): %d\n\n', n);
fprintf('Layer-wise Copper Losses:\n');
for m = 1:n
    fprintf('Layer %d: %.4f W\n', m, layer_copper_losses(m));
end

fprintf('\nTotal Copper Loss: %.4f W\n', P_cu1);
fprintf('Common Factor (I²*h/(skin_dept)*R_de): %.4f W\n\n', common_factor);

fprintf('Calculation Breakdown:\n');
fprintf('I² = %.2f A²\n', term1);
fprintf('(h/(skin_depth))*Rdc = %.6f Ω\n', term2);
fprintf('M/3 = %.3f\n', term3);
fprintf('2M² + 1 = %.1f\n\n', term4);

fprintf('Total Copper Loss (P): %.5f W\n', P_cu2);
fprintf('Equivalent to: %.2f mW\n\n', total_copper_loss2*1000);

fprintf('\nProximity Effect Ratio Analysis\n');
fprintf('P_dc: %.2f W\n', P_dc);
fprintf('Proximity Effect Ratio (F_R = P/P_dc): %.2f\n', Fr);
fprintf('Interpretation: AC losses are %.1f× greater than DC losses\n', Fr);