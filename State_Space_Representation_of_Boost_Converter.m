clc;

%-------- Model for a non-isolated boost converter

disp(["the values of plant parameters"]);
Vg = 15;
D = 0.4;
L = 2e-3;
C = 10e-6;
R = 100;

%Steady State Model of the Boost Converter
%Given by As, Bs, Cs, Ds matrices

As = [0 -(1-D)/L; (1-D)/C -1/(R*C)];
Bs = [1/L 0 0; 0 -1/C 0];
Cs = [0 1; 1 0];
Ds = [0 0 0; 0 0 0];

V0 = -Cs(1,:) * inv(As) * Bs(:,1) * Vg;
Ig = -Cs(2,:) * inv(As) * Bs(:,1) * Vg;

% Small signal model of the Boost Converter

a = [0 -(1-D)/L; (1-D)/C -1/(R*C)];
b = [1/L 0 V0/L; 0 -1/C -Ig/C];
c = [0 1; 1 0];
d = [0 0 0; 0 0 0];

ulabels = ['Vg iz d'];
ylabels = ['v0 ig'];
xlabels = ['il vc'];

printsys(As, Bs, Cs, Ds, ulabels, ylabels, xlabels);
printsys(a, b, c, d, ulabels, ylabels, xlabels);
disp(["Transfer Function in S domain"]);
disp(["Vo / d"]);
%TFB = 
TFB = tf(ss(a,b(:,3),c(1,:),[0]));
TFB
rlocus(TFB)