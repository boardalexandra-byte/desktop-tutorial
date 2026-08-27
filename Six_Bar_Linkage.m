% six bar linkage
% static equilibrium

clc;
clear;

% define the joints
A = [7 4 0];
B = [5 16 0];
C = [25 25 0];
D = [23 10 0];
E = [18 35 0];
F = [43 32 0];
G = [45 17 0];

% Define the lengths of the bars
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lBE = norm(E - B);
lEF = norm(F - E);
lFG = norm(G - F);

% Weight of Each Link
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = 0 -1 0];

% Center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3; % Not accurate
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

syms FAx FAy FBx FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx, 0, 0];
ForceC = [FCx, FCy, 0];
ForceD = [FDx, FDy, 0];
ForceE = [FEy, 0, 0];
ForceF = [FFx, FFy, 0];
ForceG = [FGx, FGy, 0];
InputTorque = [0 0 Tin];

%Applied Force
AppliedForce = [50 0 0];
% Static Equilibrium Conditions for Link AB

% Sum of Forces = 0
% Fa + Fb + Weightof AB = 0
eqn1 = ForceA + ForceB + WAB == 0;

% Sum of moments = 0
% With respect to the center of mass of link AB
% S1A x FA + S1B x FB + InputTorque = 0

eqn1 = cross(A-S1, ForceA) + cross(B-S1, ForceB) + InputTorque == 0;

% Equations for Link BEC
% Sum of Forces = 0
% -Fb + Fc + Fe + WBEC = 0
 eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;


 % Sum of Moments = 0
% With respect to the center of mass of link BEC

eqn4 = cross(B-S2, ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE)  == 0;

% -Fc + Fd + WCD = 0
% Equations for Link CD
% Sum of Forces = 0 for Link CD
eqn5 = -ForceC + ForceD + WCD == 0;

% Sum of Moments = 0 for link CD
% Sum of moments = 0 with respect to center of mass of link CD
% S3C x - FC + S3D x FD = 0
% Sum of Moments = 0 for Link CD
eqn6 = cross(C-S3, ForceC) + cross(D-S3, ForceD) == 0;

% Equations for Link EF
% Sum of Forces = 0 for Link EF
eqn7 = -ForceE + ForceF + WEF == 0;

% Sum of Moments = 0 for Link EF
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

% Equations for Link FG
% Sum of Forces = 0
eqn9 = -ForceF + ForceG + WFG + AppliedForce == 0;

% Sum of Moments = 0 for Link FG
% S5F x -FF + S5G x FG = 0
eqn10 = cross(F-S5, ForceF) + cross(G-S5, ForceG) == 0;

% Soving the 10 equations

eqnMAtric = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];

StaticSolution = solve (eqnMatric, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy]);
