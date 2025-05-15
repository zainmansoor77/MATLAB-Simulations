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
