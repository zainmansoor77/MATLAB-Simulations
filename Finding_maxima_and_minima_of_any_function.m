clc;

%{
syms x;
y = inline ('2*x^3 + 3*x^2 -12*x +17', 'x');
%grid on;
ezplot(y , [-3,2]);
grid on;
g = diff(y(x), (x));
s = solve (g)
%}

%{
clc;
syms F;
k = inline ('( F^2 * (12 - 1) ) / ( (12 * F^2)^2 + (F^2) * (F^2 - 1)^2 * (12-1)^2 * 0.5^2 )^0.5','F')
pretty(k(F))
%grid on;
%ezplot(k , F);
%grid on;
g = inline (diff(k(F), F),'F')
g
s = solve (g(F))
pretty(s)
g(-1)
%}

%{
clc;
F = linspace(0, 10, 1000);    % Define x-range
k = ( F.^2 .* (12 - 1) ) / ( (12 * F.^2).^2 + (F.^2) .* (F.^2 - 1).^2 .* (12-1).^2 .* 0.5.^2 ).^0.5;% Define the function
[max_val, idx_max] = max(k);   % Maximum value and its index
x_max = x(idx_max);            % x at which max occurs
fprintf('Maximum at x = ');
disp(x_max);
fprintf(', value = ');
disp(max_val);
%}


% Parameters
m = 6.3;
Q = 0.4;
Fx = linspace(0, 10, 1000); % Sweep Fx from 0 to 10

% Equation from image
numerator = Fx.^2 .* (m - 1);
denominator = sqrt((m .* Fx.^2 - 1).^2 + Fx.^2 .* (Fx.^2 - 1).^2 .* (m - 1).^2 .* Q^2);
K = numerator ./ denominator;

% Plotting
figure;
plot(Fx, K, 'LineWidth', 2);
grid on;
xlabel('F_x');
ylabel('K(Q, m, F_x)');
title('Plot of K(Q, m, F_x) vs F_x');