clc;
clear;

% Parameters
m = 15; % Fixed m
Q_values = 0.1:0.2:2.0; % Use fewer values for clarity
Np = 1;
Ns = 100;
Vo = 3e3;
Po_max = 300;
fr = 120e3;
Fx = linspace(0.1, 3, 10000);

% Plot setup
figure;
hold on;
grid on;
colors = lines(length(Q_values));
legend_entries = strings(length(Q_values),1);

fprintf('--- Fx values where K ≈ 1 ---\n');

for i = 1:length(Q_values)
    Q = Q_values(i);
    
    % Calculate K
    numerator = Fx.^2 .* (m - 1);
    denominator = sqrt((m .* Fx.^2 - 1).^2 + Fx.^2 .* (Fx.^2 - 1).^2 .* (m - 1).^2 .* Q^2);
    K = numerator ./ denominator;
    
    % Plot
    plot(Fx, K, 'LineWidth', 2, 'Color', colors(i,:));
    legend_entries(i) = sprintf('Q = %.1f', Q);
    
    % Find where K ≈ 1 (within 1% tolerance)
    idx_k1 = find(abs(K - 1) < 0.01, 1);
    if ~isempty(idx_k1)
        Fx_k1 = Fx(idx_k1);
        fprintf('Q = %.1f: Fx ≈ %.4f\n', Q, Fx_k1);
        plot(Fx_k1, K(idx_k1), 'ko', 'MarkerSize', 6, 'LineWidth', 1.5);
        text(Fx_k1 + 0.05, K(idx_k1), sprintf('%.2f', Fx_k1), 'FontSize', 8, 'Color', colors(i,:));
    else
        fprintf('Q = %.1f: No Fx found where K ≈ 1\n', Q);
    end
end

% Add horizontal line at K = 1
yline(1, '--k', 'K = 1', 'LabelHorizontalAlignment', 'left', 'FontWeight', 'bold');

xlabel('F_x');
ylabel('K(Q, m, F_x)');
title(sprintf('K vs F_x for different Q values at m = %.1f', m));
legend(legend_entries, 'Location', 'best');
