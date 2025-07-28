clc;

f1 = 30 * 1e3;
f2 = 20 * 1e6;
f0 = sqrt (f1 * f2);
B = f2 - f1;

Q = 0.1;

LC_ratio = 0.1;
LC = 1 / ((2 * pi * f0)^2);

%L = sqrt(LC / LC_ratio);
%L = 100 * 1e-6;
%C = L * LC_ratio;

R = 50;
L = R / (B * 2 * pi);
fprintf('L=%.4fu\n',L * 1e6);
C = 1 / ((2 * pi * f0)^2 * L);
fprintf('C=%.4fu\n',C * 1e6);
R = B * 2 * pi * L; % R calculated with the help of Bandwidth
fprintf('R=%.4f\n',R);
%R = sqrt(L / ((Q^2) * C)); % R calculated with the help of Quality factor Q
Q = 2 * pi * f0 * L / R;
fprintf('Q=%.4f\n',Q);

f1 = (1 / (2 * pi)) * [(-R / (2 *L)) + sqrt((R / (2 * L))^2 + (1 /(L * C)))];
f2 = (1 / (2 * pi)) * [(R / (2 *L)) + sqrt((R / (2 * L))^2 + (1 /(L * C)))];
fprintf('f1=%.4f\n',f1);
fprintf('f2=%.4f\n',f2);

B = f2 - f1;


fprintf('.param C=%.4fu\n',C*1e6);
fprintf('.param L=%.4fn\n',L*1e9);
fprintf('.param R=%.4f\n',R);