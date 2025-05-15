clc;

syms x;
f = inline('x^2','x');
disp(f(3));
a = diff(f(x),x);
disp(a);
%disp(class(a));

syms x;
g = inline('sin(x) / cos(x)','x');
disp(g(3));
a = diff(g(x),x);
disp(a);
%disp(class(a));
pretty(a)

syms x;
f = @(x) cos(x);
%f = inline ('cos(x)','x');
DF = diff(f(x),x);
disp(DF);
%y = inline (DF);
y = inline (DF, 'x');
disp(y(pi/2));