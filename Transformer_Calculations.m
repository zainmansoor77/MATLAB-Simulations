clc;
clear;

%{
core_loss_calculator();

function P_fe = core_loss_calculator()
    %% Input Parameters Section
    % ----------------------------
    % Material Properties
    material_type = 'Ferrite';      % Core material (Ferrite/PowderedIron)
    K_fe = 5.2e-5;                  % Core loss coefficient [W/(T^β*m³)]
    beta = 2.6;                     % Steinmetz exponent (2.6 for ferrite)
    
    % Core Geometry
    A_c = 1.2e-4;                   % Cross-sectional area [m²]
    l_m = 0.15;                     % Magnetic path length [m]
    
    % Operating Conditions
    delta_B = 0.2;                  % Peak flux density [T]
    frequency = 80e3;              % Operating frequency [Hz]
    
    %% Input Validation
    % ----------------------------
    assert(beta >= 2 && beta <= 3, 'Beta must be between 2 and 3');
    assert(A_c > 0, 'Core area must be positive');
    assert(l_m > 0, 'Path length must be positive');
    assert(delta_B > 0, 'Flux density must be positive');
    
    %% Core Volume Calculation
    % ----------------------------
    core_volume = A_c * l_m;        % [m³]
    
    %% Core Loss Calculation
    % ----------------------------
    P_fe = K_fe * (delta_B)^beta * core_volume;
    
    %% Results Presentation
    % ----------------------------
    fprintf('\nCore Loss Analysis Report\n');
    fprintf('=============================\n');
    fprintf('Material Type:       %s\n', material_type);
    fprintf('Frequency:           %.1f kHz\n', frequency/1e3);
    fprintf('Core Dimensions:\n');
    fprintf('  Cross-section:     %.2f cm²\n', A_c*1e4);
    fprintf('  Path Length:       %.1f mm\n', l_m*1e3);
    fprintf('  Volume:            %.2e m³\n', core_volume);
    fprintf('\nOperating Conditions:\n');
    fprintf('  ΔB:                %.2f T\n', delta_B);
    fprintf('  β Exponent:        %.1f\n', beta);
    fprintf('  K_fe Coefficient:  %.2e W/(T^β·m³)\n', K_fe);
    fprintf('\nCalculated Core Loss: %.2f W\n', P_fe);
    fprintf('=============================\n');
    
    %% Advanced Features (Optional)
    % ----------------------------
    % Uncomment to generate warning at 200mT saturation limit
    % if delta_B > 0.2
    %     warning('Flux density exceeds typical ferrite saturation limit');
    % end
end

%}


%% Transformer Design Calculator (Procedural Version)
% Implements Chapter 12 equations from "Fundamentals of Power Electronics"

function transformer_design()
    %% Example Parameters 
    params.core_material = 'Ferrite';
    params.Ac_cm2 = 0.635;       % 2213 pot core
    params.WA_cm2 = 0.297;
    params.MLT_cm = 4.5;
    params.lm_cm = 3.2;
    params.B_sat = 0.35;
    params.Kfe = 24.7;
    params.beta = 2.6;
    params.rho = 1.724e-6;
    params.Ku = 0.5;
    params.lambda1 = 62.5e-6;
    params.Itot = 8;
    params.Ptot_max = 0.25;
    params.turns_ratios = 1/5;

    %% Run Design Calculations
    [design, status] = optimize_design(params);
    
    %% Display Results
    if status == 1
        print_report(design);
    else
        fprintf('Design failed: %s\n', design.error);
    end
end

function [design, status] = optimize_design(params)
    status = 1; % Assume success
    design = params;
    
    try
        %% Calculate Core Constants (Eq 12.16)
        term1 = (design.WA_cm2 * design.Ac_cm2^(2*(design.beta-1)/design.beta)) / ...
               (design.MLT_cm * design.lm_cm^(2/design.beta));
        term2 = (design.beta/2)^(-design.beta/(design.beta+2)) + ...
                (design.beta/2)^(2/(design.beta+2));
        design.Kgfe = term1 * term2;

        %% Verify Core Size (Eq 12.19)
        Kgfe_req = (design.rho * design.lambda1^2 * design.Itot^2 * ...
                   design.Kfe^(2/design.beta)) / ...
                   (4 * design.Ku * design.Ptot_max^((design.beta+2)/design.beta)) * 1e8;
        
        if design.Kgfe < Kgfe_req
            warning('Core is undersized - select larger core');
        end

        %% Optimal Flux Density (Eq 12.20)
        numerator = 1e8 * design.rho * design.lambda1^2 * design.Itot^2 * design.MLT_cm;
        denominator = 2 * design.Ku * design.WA_cm2 * design.Ac_cm2^3 * design.lm_cm * ...
                     design.beta * design.Kfe;
        design.deltaB = (numerator / denominator)^(1/(design.beta+2));

        %% Safety Checks
        if design.deltaB > 0.75 * design.B_sat
            error('Flux density exceeds 75%% of saturation limit');
        end

        %% Primary Turns (Eq 12.21)
        design.n1 = (design.lambda1 * 1e4) / (2 * design.deltaB * design.Ac_cm2);

        %% Winding Calculations
        design.n_secondary = design.n1 * design.turns_ratios;
        design.alpha = (design.n_secondary .* design.turns_ratios) / ...
                      (design.n1 * design.Itot);
        design.wire_areas = (design.alpha * design.Ku * design.WA_cm2) ./ ...
                           [design.n1, design.n_secondary];

        %% Loss Calculations
        design.Pfe = design.Kfe * design.deltaB^design.beta * ...
                    (design.Ac_cm2 * design.lm_cm);
        Pcu_factor = (design.rho * design.lambda1^2 * design.Itot^2) / ...
                    (4 * design.Ku) * (design.MLT_cm / (design.WA_cm2 * design.Ac_cm2^2));
        design.Pcu = Pcu_factor / design.deltaB^2;

    catch ME
        design.error = ME.message;
        status = 0;
    end
end

function print_report(design)
    fprintf('\n=== Transformer Design Report ===\n');
    fprintf('Core Material:          %s\n', design.core_material);
    fprintf('Core Size:              %.2f cm²\n', design.Ac_cm2);
    fprintf('Window Area:            %.2f cm²\n', design.WA_cm2);
    fprintf('\n--- Optimal Operating Point ---\n');
    fprintf('ΔB:                     %.3f T\n', design.deltaB);
    fprintf('Primary Turns (n1):     %.1f\n', design.n1);
    fprintf('Secondary Turns:        %s\n', mat2str(design.n_secondary,2));
    fprintf('\n--- Loss Analysis ---\n');
    fprintf('Core Loss:              %.2f W\n', design.Pfe);
    fprintf('Copper Loss:            %.2f W\n', design.Pcu);
    fprintf('Total Loss:             %.2f W\n', design.Pfe + design.Pcu);
    fprintf('Allowed Total Loss:     %.2f W\n', design.Ptot_max);
    fprintf('\n--- Safety Margins ---\n');
    fprintf('Saturation Margin:      %.1f%%\n', ...
           (design.B_sat/design.deltaB - 1)*100);
    fprintf('=================================\n');
end



transformer_design();