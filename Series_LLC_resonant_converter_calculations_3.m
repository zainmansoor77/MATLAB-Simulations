clc;

% Parameters
lemda = 0.1667; % m = 6

%Q = 0.105; % previous val = 0.6

Np = 1;           % Primary turns
Ns = 75;          % Secondary turns
n = Np / Ns;
Vo = 15*1e3;      % Output voltage in Volts
Po_max = 300;     % Max output power in Watts
%Qmax = Q;        % Quality factor (already given)
fr = 10e3;        % Resonant frequency (Hz)
fsw = 50e3;       % Switching Frequency (Hz)
R0 = 30e3;
Cr = 4.7 * 1e-6;

Lr = ((1 / (2 * pi * fr)).^2) / Cr;
% Display results
fprintf('Selected Cr = %.3e F\n', Cr);
fprintf('Calculated Lr primary = %.3e H\n', Lr);

freq = (1 / (2 * pi * sqrt(Lr * Cr)) );
fprintf('Calculated Fr = %.1f KHz\n', freq * 1e-3);

Z0 = sqrt(Lr / Cr); % characteristic impedance 
fprintf('Calculated Z0 = %.3f Ohms\n', Z0);

% Eq. 4: Calculate Rac_min
%Rac_min = (8/pi^2) * (Np^2 / Ns^2) * (Vo^2 / Po_max);
Rac_min = (8/pi^2) * (Np^2 / Ns^2) * R0;
fprintf('Rac_min = %.2f Ohms\n', Rac_min);

Q = Z0 / Rac_min;

fprintf('Qmax = %.4f\n',Q);

