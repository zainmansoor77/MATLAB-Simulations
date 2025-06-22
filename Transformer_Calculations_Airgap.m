clc;

Turn_ratio = 75;                       % Ns / Np = 75
Lm = 526.5557 * 1e-6;                  % Lm = 526.5557uH
lamda = 1.81913  * 1e-3;               % lamda = 1.81913mV-s
delta_B = 0.25;                         % 0.2T peak-to-peak flux density
min_Np = 0;                            % Minimum turns to avoid core saturation
Np = 0;                                % Number of turns in primary side
Lm_ungapped = 0;                       % Ungapped inductance at the Minimum turns, to avoid core saturation
lm_ungapped_Le = 0;
%Lm_ungapped2 = 0;                     % Ungapped inductance at greater primary number of turns
count = 0;                             % By default it should be 0, otherwise it would skip the part that verify that no air gap is required
U0 = 4 * pi * 1e-7;                    % Permeability of free space [H/m]
Lg = 0.00;

% %{

% EE160/80/40 core parameters (N87 material)

Ae = 1.64 * 1e-3;                      % Effective cross-sectional area = 1640mm^2 or (1.64 * 1e-3)m^2
le = 0.381;                            % Effective magnetic length path = 381mm or 0.381m
Ve = 6.24 * 1e-4;                      % Effective core volume = 623930mm^3 or (6.24 * 1e-4)m^3
AL = 10700 * 1e-9;                     % Ungapped AL(N87 material) = 10700nH / turn^2
Fsw = 50 * 1e3;                        % Switching Frequency of 50KHz
Duty_Cycle = 0.5;                      % Assumed full-bridge square wave (duty = 50%).
%Ur = 1970;


% %}

%{

% U93/76/30 core parameters (3C90 or 3C94 material)

Ae = 0.84 * 1e-3;                      % Effective cross-sectional area = 840mm^2 or (0.84 * 1e-3)m^2
le = 0.354;                            % Effective magnetic length path = 354mm or 0.354m
Ve = 2.97 * 1e-4;                      % Effective core volume = 297000mm^3 or (2.97 * 1e-4)m^3
AL = 6400 * 1e-9;                      % Ungapped AL(N87 material) = 6400nH / turn^2
Fsw = 50 * 1e3;                        % Switching Frequency of 50KHz
Duty_Cycle = 0.5;                      % Assumed full-bridge square wave (duty = 50%).
%Ur = 1970;


%}

%{

% U126/91/20 core parameters (3C90 or 3C94 material)

Ae = 0.56 * 1e-3;                     % Effective cross-sectional area = 560mm^2 or (0.56 * 1e-3)m^2
le = 0.48;                            % Effective magnetic length path = 480mm or 0.48m
Ve = 2.688 * 1e-4;                    % Effective core volume = 268800mm^3 or (2.688 * 1e-4)m^3
AL = 3000 * 1e-9;                     % Ungapped AL(N87 material) = 3000nH / turn^2
Fsw = 50 * 1e3;                       % Switching Frequency of 50KHz
Duty_Cycle = 0.5;                     % Assumed full-bridge square wave (duty = 50%).
%Ur = 1970;


%}

%% Minimum number of turns to avoid core saturation

min_Np = lamda / (delta_B * Ae);
fprintf('Minimum number of Primary turns required: %.2f \n', min_Np);

%% Verify Magnetizing Inductance

Lm_ungapped = AL * (min_Np ^ 2);
fprintf('Magneuizing inductance at the Minimum Np: %.3f uH \n', Lm_ungapped * 1e6);

Ur = (AL * le) / (U0 * Ae);
fprintf('Ur for this material is %d \n', Ur);
lm_ungapped_Le = (U0 * Ur * (min_Np ^ 2) * Ae) / le;
fprintf('Magneuizing inductance at the Minimum Np using Le: %.3f uH \n', lm_ungapped_Le * 1e6);

fprintf('\n');

for i = ceil(min_Np):100
    if (Lm - Lm_ungapped) < 0
        break;
    end
    
    Lm_ungapped = AL * (i ^ 2);
    fprintf('Np = %d ... Magnetizing inductance at this Np: %.3f uH \n', i, Lm_ungapped * 1e6);

end

fprintf('\n');

if ( abs(Lm - (AL * ((i - 2) ^ 2))) < (20 * 1e-6) ) || ( abs(Lm_ungapped - Lm) < (20 * 1e-6) )
    fprintf('No Air Gap required \n');
else
    count = 1;
    fprintf('Air Gap required \n');
end

fprintf('\n');

if (count == 0) && ( (Lm - (AL * ((i - 2) ^ 2))) < (20 * 1e-6) )
    Np = i - 2;
    fprintf('Np = %d ... Magnetizing inductance at this Np: %.3f uH \n', Np, (AL * ((i - 2) ^ 2)) * 1e6);

elseif (count == 0) && ( (Lm_ungapped - Lm) < (20 * 1e-6) )
    Np = i - 1;
    fprintf('Np = %d ... Magnetzing inductance at this Np: %.3f uH \n', Np, Lm_ungapped * 1e6);
end

if count == 1
    fprintf('\n');
    Np = i -1;
    fprintf('Np = %d ... Magnetizing Ungapped inductance at this Np: %.3f uH \n', Np, Lm_ungapped * 1e6);
    delta_B = lamda / (Np * Ae);
    fprintf('The Maximum B we will get is: %.2f T \n', delta_B);
    %Lg = (((Np ^ 2) * U0 * Ae) / Lm);
    Lg = (((Np ^ 2) * U0 * Ae) / (Lm)) - (le / Ur);
    fprintf('The air gap required is %.03f mm \n', Lg * 1e3);
    
    fprintf('\n');
end