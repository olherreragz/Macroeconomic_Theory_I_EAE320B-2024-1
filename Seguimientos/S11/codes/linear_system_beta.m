clear all;

a1 = 0.6*0.96;
a2 = 0.15*0.96*(0.9/0.97)^2;
a3 = 0.15*0.96*(0.9/1)^2;
a4 = 0.1*0.96*(0.9/1.03)^2;
b1 = 0,9;
b2 = 0.97;
b3 = 1;
b4 = 1.03;

c1 = 0.05*0.96*(0.97/0.9)^2;
c2 = 0.65*0.96;
c3 = 0.25*0.96*(0.97)^2;
c4 = 0.05*0.96*(0.97/1.03)^2;
d1 = 0,9;
d2 = 0.97;
d3 = 1;
d4 = 1.03;

e1 = 0.01*0.96*(1/0.9)^2;
e2 = 0.07*0.96*(1/0.97)^2;
e3 = 0.085*0.96;
e4 = 0.07*0.96*(1.03/0.97)^2;
f1 = 0,9;
f2 = 0.97;
f3 = 1;
f4 = 1.03;

g1 = 0.015*0.96*(1.03/0.9)^2;
g2 = 0.05*0.96*(1.03/0.97)^2;
g3 = 0.28*0.96*(1.03)^2;
g4 = 0.655*0.96;
h1 = 0,9;
h2 = 0.97;
h3 = 1;
h4 = 1.03;

syms p1 p2 p3 p4
eqn1 = a1*(p1 + b1) + a2*(p1 + b2) + a3*(p1 + b3) + a4*(p1 + b4) == p1;
eqn2 = c1*(p2 + d1) + c2*(p2 + d2) + c3*(p2 + d3) + c4*(p2 + d4) == p2;
eqn3 = e1*(p3 + f1) + e2*(p3 + f2) + e3*(p3 + f3) + e4*(p3 + f4) == p3;
eqn4 = g1*(p4 + h1) + g2*(p4 + h2) + g3*(p4 + h3) + g4*(p4 + h4) == p4;

sol = solve([eqn1, eqn2, eqn3, eqn4], [p1, p2, p3, p4]);
p1 = double(sol.p1)
p2 = double(sol.p2)
p3 = double(sol.p3)
p4 = double(sol.p4)


% [A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [p1, p2, p3, p4]);
% X = linsolve(A,B);


% Define the Markov matrix
P = [0.6 0.15 0.15 0.1;
     0.05 0.65 0.25 0.05;
     0.01 0.07 0.85 0.07;
     0.015 0.05 0.28 0.655];

% Define the symbolic variables for the stationary distribution
syms pi1 pi2 pi3 pi4

% Set up the system of linear equations
eqn1 = pi1*P(1,1) + pi2*P(2,1) + pi3*P(3,1) + pi4*P(4,1) == pi1;
eqn2 = pi1*P(1,2) + pi2*P(2,2) + pi3*P(3,2) + pi4*P(4,2) == pi2;
eqn3 = pi1*P(1,3) + pi2*P(2,3) + pi3*P(3,3) + pi4*P(4,3) == pi3;
eqn4 = pi1*P(1,4) + pi2*P(2,4) + pi3*P(3,4) + pi4*P(4,4) == pi4;
eqn5 = pi1 + pi2 + pi3 + pi4 == 1;

% Solve the system of equations
sol = solve([eqn1, eqn2, eqn3, eqn4, eqn5], [pi1, pi2, pi3, pi4]);

% Display the results
pi1 = double(sol.pi1)
pi2 = double(sol.pi2)
pi3 = double(sol.pi3)
pi4 = double(sol.pi4)

Expected_value = pi1 * p1 + pi2 * p2 + pi3 * p3 + pi4 * p4;


