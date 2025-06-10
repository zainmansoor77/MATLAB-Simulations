clc;

K_fe = 24.7;     % Kf e is a constant of proportionality which depends on the operating frequency. (Core loss coefficient) ( material constant)
K_gfe = 30.4 * 1e-3; % It is the core geometrical constant Kgfe (cm^x)
Kgfe_req = 0; % minimum required K_gfe (cm^x)
delta_B = 1;  % peak ac flux density (ΔB) in Tesla
beta = 2.6;     % material constant
%Ac = 1;       % cross-sectional core area (Ac) in cm²
lm = 10.3;       % core length (lm) in cm
no_w = 2;     % number of windings in a transformer
n1 = 2.10;       % Primary turns n1
n2 = 210;       % Secondary turns n2
n = [n1, n2]; % number of turns in each winding
turn_ratio = [100]; % turns ratio (Ns / Np)
lamda1 = 1.79213 * 1e-3;   % Volts second blance in V-sec
I_rms1 = 1.9854;   % RMS current of the primary winding
I_rms2 = 19.494 * 1e-3;   % RMS current of the secondary winding
I_rms = [I_rms1, I_rms2]; % RMS current of each winding
I_tot = 0;    % It is the sum of the rms winding currents, referred to winding 1
I_ref_pri = zeros (1, no_w);  % It is the rms winding currents of each winding, referred to winding 1
Imag_peak = 0;% Peak ac magnetizing current, referred to winding 1:
rho = 1.724 * 1e-9;      % resistivity of the winding coil
Ku = 0.8;       % Winding fill factor
MLT = 7.62;      % Mean Length per Turn [cm]
Wa = 2.13;       % Core winding area cm²
Ac = 1.74;       % Ac is the core cross-sectional area  cm²
P_fe = 0;     % Core loss of the transformer
P_cu = 0;     % Copper loss of the transformer
P_tot1 = 10;   % Total Power loss in W (P_cu + P_fe)
P_tot2 = 10;   % Total power loss in W (optimal total power loss by the complex equation)
Alpha = zeros (1, no_w);    % fraction of the window area allocated to each winding
Aw = [52.41*1e-3, 3.243*1e-3];      % wire area (guage) cm²
Aw_req = zeros (1, no_w);   % Wire area minimum requirement cm²
L_mag = 1.89 * 1e-3;   % magnetizing inductance, referred to winding 1
Ll1 = 172 * 1e-6;
Ll2 = 1.72;
U0 = 4 * pi * 1e-7; % Permeability of free space [H/m]
Ur = 2000;       % relative permeability
U = U0 * Ur;        % Core permeability [H/m]
Rw = zeros (1, no_w); % Winding resistance

L12_choice = 2; % 1 for  L12 = (n1 * n2) / R; 2 for L12 = (n2 / n1) * L_M;

% Calculate I total
for m = 1:no_w
    I_ref_pri(m) = (n(m) / n(1)) * I_rms(m);
end

I_tot = sum(I_ref_pri);

fprintf('Total sum of the rms winding currents, referred to winding 1: %.4f A\n', I_tot);

% Calculate peak ac flux density (ΔB) in Tesla
term1 = (rho  * lamda1^2 * I_tot^2) / (2 * Ku); % specification regarding the applications (ρ, Itot, λ1, Ku, Ptot)
term2 = (MLT) / (Wa * (Ac^3) * lm);             % function of the core geometry
term3 = 1 / (beta * K_fe);
delta_B = (1e8 * term1 * term2 * term3)^(1 / (beta + 2));

fprintf('\npeak ac flux density (ΔB): %.4f T\n', delta_B);
fprintf('<strong>Check for core saturation </strong>\n');


% Calculate primary turns n1
%n1 = (lamda1 / (2 * delta_B * Ac)) * 1e4;
n1 = sqrt((L_mag * lm) / (U * Ac));
n(1) = n1;
fprintf('primary turns n1: %.4f \n', n(1));

for m = 2:no_w 
    n(m) = turn_ratio(m-1) * n(1) ;
    fprintf('secondary turns n%.0f: %.4f \n',m, n(m));
end



% Calculate core loss
P_fe = K_fe * (delta_B)^beta * Ac * lm; % Eq. 12.1

fprintf('Total Core Loss (P_fe): %.4f W\n', P_fe);

% Calculate copper loss
term4 = (rho  * lamda1^2 * I_tot^2) / (4 * Ku); % specification regarding the applications (ρ, Itot, λ1, Ku, Ptot)
term5 = (MLT) / (Wa * Ac^2);                    % function of the core geometry
term6 = (1 / delta_B)^2;                        % function of delta_B

P_cu = term4 * term5 * term6 * 1e8;
%fprintf('specification regarding the applications: %.4f p\n', term4*1e12);
%fprintf('function of the core geometry: %.4f u\n', term5 * 1e6);
%fprintf('function of delta_B: %.4f u\n', term6 * 1e6);
fprintf('Total Copper Loss (P_cu): %.4f W\n', P_cu);

% Calculate total power (copper + core) loss
P_tot1 = P_cu + P_fe;
fprintf('Total Power (copper + core) Loss (P_total): %.4f W\n', P_tot1);

% Calculate total power (from complex equation) loss
term7 = (Ac * lm * K_fe)^(2 / (beta + 2));
term8 = (((rho * (lamda1^2) * (I_tot^2)) * MLT) / (4 * Ku * Wa * (Ac^2)))^(beta / (beta + 2));
term9 = ((beta / 2)^(-(beta / (beta + 2)))) + ((beta / 2)^(2 / (beta + 2)));
P_tot2 = term7 * term8 * term9 * 1e8;
fprintf('Total Power (from complex equation) Loss (P_total): %.4f W\n', P_tot2);

% Calculate Kgfe (Core geometrical constant)
%term10 = (Wa * (Ac)^(2*(beta-1)/beta)) / (MLT * lm^(2/beta));
%term11 = ((beta/2)^(-beta/(beta+2)) + (beta/2)^(2/(beta+2)))^(-(beta + 2) / beta);
%K_gfe = (term10 * term11) * 1e8;
%fprintf('The K_gfe constant: %.4f \n', K_gfe);

%calculate the K_gfe requirement
Kgfe_req = ((rho * (lamda1^2) * I_tot^2 * K_fe^(2/beta)) / (4 * Ku * P_tot1^((beta+2)/beta))) * 1e8;
fprintf('The K_gfe requirement constant: %.4f \n', Kgfe_req);
        
if K_gfe < Kgfe_req
    fprintf('<strong>Core is undersized - select larger core </strong>\n');
else
    fprintf('<strong>Core is of correct size </strong>\n');
end

%fprintf('Hello <strong> bold </strong> world.\n');  % this is how we print statement in bold in matlab script  

%{
Kgfe is a measure of the magnetic size of a core, for applications in which core loss is significant. Unfortunately, Kgfe depends on β, and hence
the choice of core material affects the value of Kgfe. However, the β of most high-frequency ferrite materials lies in the narrow range 2.6 to 2.8, 
and Kgfe varies by no more than ±5% over this range. Appendix B lists the values of Kgfe for various standard ferrite cores, for the value β = 2.7.

Once a core has been selected, then the values of Ac, WA, m, and MLT are known. The peak ac flux density ΔB can then be evaluated, , and the primary
turns n1 can be found. The number of turns for the remaining windings can be computed using the desired turns ratios. The various window area 
allocations are found.
%}

% Calculate the fraction of window area allocated to each winding.
for m = 1:no_w
    Alpha(m) = (n(m) * I_rms(m)) / (n(1) * I_tot);
    fprintf('fraction window area allocated for n%.0f winding: %.4f \n', m, Alpha(m));

    Aw_req(m) = (Alpha(m) * Ku * Wa) / n(m);
    fprintf('Wire size requirement for n%.0f winding: %.4f \n',m, Aw_req(m));

    if Aw(m) <= Aw_req(m)
        fprintf('<strong>Wire guage is fine for n%.0f</strong>\n', m); 
    else
        fprintf('<strong>Reduce the AWG for n%.0f</strong>\n', m); 
    end
end

U = Ur * U0; % Core permeability [H/m]
L_mag = (U * (n1^2) * Ac) / lm;
fprintf('Magnetizing inductance, referred to winding 1: %.4f H (%.3f mH) \n', L_mag, L_mag * 1e3);

Imag_peak = lamda1 / (2 * L_mag);
fprintf('Peak ac magnetizing current, referred to winding 1: %.4f A (%.3f uA) \n', Imag_peak, Imag_peak * 1e6);

for m = 1:no_w
    Rw(m) = (rho * n(m) * MLT) / Aw(m);
    fprintf('Winding Resistance of the n%.0f winding: %.4f Ohms \n', m, Rw(m));
end


%% Calculations for Real Transformer

% 1. Compute core reluctance (Eq. 10.34)
R = (lm * 1e-2) / (U * (Ac * 1e-4)); 
fprintf('Core Reluctance (R): %.4e A-turns/Wb\n', R);

% 2. Compute total MMF (Eq. 10.35)
MMF = n1 * I_rms1 + n2 * I_rms2;
fprintf('Total MMF (n1i1 + n2i2): %.2f A-turns\n', MMF);
    
% 3. Compute flux (Eq. 10.36)
Phi = MMF / R;
fprintf('Magnetic Flux (Φ): %.4e Wb\n', Phi);

% 2. Compute magnetizing inductance (Eq. 10.44)
L_M = (n1^2) / R;
fprintf('Magnetizing Inductance calculated for 2nd time (L_M): %.4f mH\n', L_M * 1e3);
    
% 3. Compute magnetizing current (Eq. 10.44)
i_M = I_rms1 + (n2/n1) * I_rms2;
fprintf('Magnetizing current (i_M): %.4f A\n', i_M);

if L12_choice == 1
    L12 = (n1 * n2) / R; % Eq. 10.49 (using reluctance)
elseif L12_choice == 2
     L12 = (n2 / n1) * L_mag; % Eq. 10.49 (using L_M)
end



% Self-inductances (Eq. 10.50)
L11 = Ll1 + (n1 / n2) * L12;
L22 = Ll2 + (n2 / n1) * L12;

% Effective turns ratio (Eq. 10.51)
n_e = sqrt(L22 / L11);

% Coupling coefficient (Eq. 10.52)
k = L12 / sqrt(L11 * L22);

fprintf('Mutual Inductance (L12): %.9f H (%.3f uH)\n', L12, L12 * 1e6);
fprintf('Primary Self-Inductance (L11): %.9f H (%.3f uH)\n', L11, L11 * 1e6);
fprintf('Secondary Self-Inductance (L22): %.9f H (%.3f uH)\n', L22, L22 * 1e6);
fprintf('Effective Turns Ratio (n_e): %.4f\n', n_e);
fprintf('Coupling Coefficient (k): %.4f\n', k);