clc;

% Parameters
lambda = 0.2; % m = 8

Q = 0.1; % previous val = 0.6

Np = 1;           % Primary turns
Ns = 75;          % Secondary turns
n = Np / Ns;
Vo = 15*1e3;      % Output voltage in Volts
Po_max = 500;     % Max output power in Watts
%Qmax = Q;        % Quality factor (already given)
fr = 20e3;        % Resonant frequency (Hz)
fsw = 50e3;       % Switching Frequency (Hz)
R0 = 30e3;
Cr = 0.5 * 1e-6;
M_max = 1.0;
M_min = 0.75;

% Eq. 4: Calculate Rac_min
Rac_min = (8/pi^2) * (Np^2 / Ns^2) * (Vo^2 / Po_max);
%Rac_min = (8/pi^2) * (Np^2 / Ns^2) * R0;
fprintf('Rac_min = %.2f Ohms\n', Rac_min);

%Lr = ((1 / (2 * pi * fr)).^2) / Cr;
Lr = ((Qmax * Rac_min).^2) * Cr; 
% Display results
fprintf('Selected Cr = %.3e F\n', Cr);
fprintf('Calculated Lr primary = %.3e H\n', Lr);

freq = (1 / (2 * pi * sqrt(Lr * Cr)) );
fprintf('Calculated Fr = %.1f KHz\n', freq * 1e-3);

Z0 = sqrt(Lr / Cr); % characteristic impedance 
fprintf('Calculated Z0 = %.3f Ohms\n', Z0);


%Q = Z0 / Rac_min;
fprintf('Qmax = %.4f\n',Q);

Lm = Lr / lambda;
fprintf('Lm primary = %.2e H\n', Lm);

% Calculate turns ratio
a = Np / Ns;
% Calculate secondary side inductances
Lm_secondary = Lm / (a^2);
Llk_secondary = Lr / (a^2);
% Display results
fprintf('Secondary Magnetizing Inductance: %.3e H\n', Lm_secondary);
fprintf('Secondary Leakage Inductance:     %.3e H\n', Llk_secondary);

fn = fsw / fr; % Normalized frequency
fprintf('Calculated Fn = %.1f\n', fn);

% Calculate voltage gain
term1 = (1 + lambda - (lambda./(fn^2))).^2;
term2 = (Q^2) * (fn - 1./fn).^2;
M = 1 ./ sqrt(term1 + term2);
fprintf('Gain of the Resonant tank: %.3f \n', M);

% Calculate the no load gain
Mol = 1 / (1 + lambda - (lambda / (fn^2)));
fprintf('Open Load Gain of the Resonant tank: %.3f \n', Mol);

term3 = 1 / (M_max^2);
term4 = (1 + lambda - (lambda./(fn^2))).^2;
term5 = (fn - 1./fn).^2;
Q_min = sqrt((term3 - term4) / term5);
fprintf('Q_min : %.3f \n', Q_min);

term3 = 1 / (M_min^2);
term4 = (1 + lambda - (lambda./(fn^2))).^2;
term5 = (fn - 1./fn).^2;
Q_max = sqrt((term3 - term4) / term5);
fprintf('Q_max : %.3f \n', Q_max);



%% Plotting

Q_values = [0,0.1, Q, 0.5, 1, 2, 5];  % Quality factors to analyze
fn_min = 0.2;
fn_max = 6.0;
fn = linspace(fn_min, fn_max, 10000); % Create frequency vector
M = zeros(length(Q_values), length(fn)); % Initialize gain matrix

% Calculate voltage gain for each Q
for i = 1:length(Q_values)
    Q2 = Q_values(i);
    
    % Calculate denominator components
    term1 = (1 + lambda - (lambda./(fn.^2))).^2;
    term2 = (Q2^2) * (fn - 1./fn).^2;
    
    % Calculate voltage gain
    M(i,:) = 1 ./ sqrt(term1 + term2);
end

% Plot results
figure('Position', [100, 100, 800, 600], 'Color', 'white');
% Create colormap for different Q values
colors = lines(length(Q_values));
hold on;
grid on;

% Plot gain curves
for i = 1:length(Q_values)
    plot(fn, M(i,:), 'LineWidth', 2.5, ...
        'Color', colors(i,:), ...
        'DisplayName', sprintf('Q2 = %.1f', Q_values(i)));
end

% Add special markers and lines
plot([1, 1], [0, max(M(:))], 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Series Resonance (f_n = 1)');

% Calculate parallel resonance frequency (f_no)
f_no = 1/sqrt(1 + 1/lambda);
plot([f_no, f_no], [0, max(M(:))], 'm--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Parallel Resonance (f_n = %.2f)', f_no));

% Label peak gain points
for i = 1:length(Q_values)
    [M_max, idx] = max(M(i,:));
    f_peak = fn(idx);
    plot(f_peak, M_max, 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colors(i,:), ...
        'HandleVisibility', 'off');
    
    % Add text label
    text(f_peak, M_max+0.05, sprintf('(%.2f, %.2f)', f_peak, M_max), ...
        'FontSize', 9, 'Color', colors(i,:), ...
        'HorizontalAlignment', 'center');
end

% Plot formatting
title(sprintf('LLC Resonant Converter Voltage Gain (\\lambda = %.1f)', lambda), ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel('Normalized Frequency (f_n)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Voltage Gain (M)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
xlim([fn_min, fn_max]);
ylim([0, min(4, max(M(:)) + 0.2)]);

% Set logarithmic scale for better visualization
set(gca, 'XScale', 'log', 'FontSize', 11);

% Add grid
grid minor;
box on;

% Key Observations Annotation
annotation('textbox', [0.15, 0.65, 0.2, 0.2], 'String', {
    'Key Observations:',
    sprintf('• Series resonance at f_n = 1 (M = 1)'),
    sprintf('• Parallel resonance at f_n = %.2f', f_no),
    '• Gain peaks below resonance',
    '• Higher Q = narrower bandwidth',
    '• Lower Q = higher peak gain'
    }, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', [1 1 1 0.7], ...
    'EdgeColor', 'none', ...
    'FontSize', 10);

%% Calculation

fprintf('.param Cr=%.4fu\n',Cr*1e6);
fprintf('.param Lr=%.4fu\n',Lr*1e6);
fprintf('.param Lm=%.4fu\n',Lm*1e6);
fprintf('.param LSm=%.4fm\n',Lm_secondary*1e3);
fprintf('.param LSr=%.4fm\n',Llk_secondary*1e3);