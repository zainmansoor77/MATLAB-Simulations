clc;

%equation = [1 4];
%roots(equation)

syms x
equation = x + 4 == 0;
sol = solve([equation],[x]);
fprintf('sol = %d \n',sol);
%disp(sol);

equation = [1 0 4];
roots(equation)

equation = [1 0 -4];
roots(equation)

equation = [2 1 0 4];
roots(equation)

syms x y;
eq1 = 2*x + 8*y == 15;
eq2 = 9*x - 6*y == 21;
sol = solve([eq1,eq2],[x,y]);
fprintf('sol x = ')
disp(sol.x);
fprintf('\nsol y = ')
disp(sol.y);

syms x y z;
eq1 = 2*x + 8*y - 9*z == 15;
eq2 = 9*x - 6*y + 5*z == 21;
eq3 = 12*x + 19*y - 13*z == 26;
sol = solve([eq1,eq2,eq3],[x,y,z]);
fprintf('sol x = ')
disp(sol.x);
fprintf('\nsol y = ')
disp(sol.y);
fprintf('\nsol z = ')
disp(sol.z);