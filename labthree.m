clc
clear
close all
% Mechanical bathroom scale - lumped model referred to the pinion angle theta

%%  Unit conversions 
in2m = 0.0254;
mm2m = 1e-3;
g    = 9.81;

%% Load 
M_person = 70;              % kg (weight being measured)
W        = M_person*g;      % N

%% Lengths (Table 1)
L_CD = 5.2  * in2m;         % short lever
L_AD = 4.85 * in2m;         % long lever pivot A -> D
L_AB = 9.2  * in2m;         % long lever pivot A -> nose B
L_BP = 0.25 * in2m;         % lever arm input arm
L_PS = 0.65 * in2m;         % lever arm output arm (to rack)
L_AE = 2.42 * in2m;         % long lever pivot A -> platform knife edge E
L_DF = 2.3  * in2m;

%% Inertias and masses (Table 1) 
J_CD   = 0.000176;          % each short lever, about C
J_AB   = 0.0002;            % each long lever, about A
J_PB   = 0.000123;          % lever arm parts, about P
J_PS   = 0.000123;
M_B    = 0.01;              % nose bracket at B
M_rack = 0.046754;

% pinion and dial
R_pinion = 0.075 * in2m; 
t_pinion = 0.096 * in2m; 
rho_pin  = 7850;            % steel (2700 for aluminum)
m_pinion = rho_pin*pi*R_pinion^2*t_pinion;
J_pinion = 0.5*m_pinion*R_pinion^2;


m_dial = 0.015;             % kg
r_dial = 0.05;              % m
J_dial = 0.5*m_dial*r_dial^2;

% springs
G = 81.7e9;  E = 203.4e9;

k1 = G*(1.3*mm2m)^4 / (8*(6.2*mm2m)^3*9);   % each main spring, N/m
K1 = 2*k1;                                   % two in parallel

K_B = E*(2.0*mm2m)^4 / (67.8*(13*mm2m)*6);   % torsion spring, N*m/rad

mT  = 0.44e-3;
K_T = G*(0.40*mm2m)^4 / (8*(2.4*mm2m)^3*40); % rack spring, N/m

%kinematics 
G_L  = L_PS/L_BP; 
r_S  = R_pinion; 
r_B  = R_pinion/G_L; 
r_E  = r_B*L_AE/L_AB;
phiP  = R_pinion/L_PS;          
phiAB = r_B/L_AB;              
phiCD = phiAB*L_AD/L_CD;      

%parameters
J_eq = J_pinion + J_dial ...
     + (M_rack + mT/3)*r_S^2 ...
     + (J_PB + J_PS)*phiP^2 ...
     + 2*J_AB*phiAB^2 + 2*J_CD*phiCD^2 ...
     + M_B*r_B^2 ...
     + M_person*r_E^2;             % person moves with the platform

K_eq = K1*r_B^2 + K_B*phiP^2 + K_T*r_S^2;

T_in = W*r_E; 

% damping
zeta = 0.15; 
B_eq = 2*zeta*sqrt(K_eq*J_eq);

% differential equation
syms theta(t)
eqn1   = J_eq*diff(theta,t,2) + B_eq*diff(theta,t) + K_eq*theta == T_in;
Dtheta = diff(theta,t);
ICs    = [theta(0) == 0, Dtheta(0) == 0];

solutionTheta = dsolve(eqn1, ICs);
solutionOmega = diff(solutionTheta, t);
solutionAlpha = diff(solutionOmega, t);
solutionRead  = solutionTheta*K_eq/(r_E*g);   % dial reading in kg

tEnd = 3;

%Plots
figure
fplot(solutionTheta, [0 tEnd], 'LineWidth', 1.5, 'MeshDensity', 2000)
grid on
title('Pinion angle vs time', 'Interpreter', 'latex')
xlabel('Time, t (s)', 'Interpreter', 'latex')
ylabel('Pinion angle, $\theta(t)$ (rad)', 'Interpreter', 'latex')

figure
fplot(solutionOmega, [0 tEnd], 'LineWidth', 1.5, 'MeshDensity', 2000)
grid on
title('Pinion angular velocity vs time', 'Interpreter', 'latex')
xlabel('Time, t (s)', 'Interpreter', 'latex')
ylabel('Angular velocity, $\omega(t) = \dot{\theta}(t)$ (rad/s)', 'Interpreter', 'latex')

figure
fplot(solutionAlpha, [0 tEnd], 'LineWidth', 1.5, 'MeshDensity', 2000)
grid on
title('Pinion angular acceleration vs time', 'Interpreter', 'latex')
xlabel('Time, t (s)', 'Interpreter', 'latex')
ylabel('Angular acceleration, $\alpha(t) = \ddot{\theta}(t)$ (rad/s$^2$)', 'Interpreter', 'latex')

figure
fplot(solutionRead, [0 tEnd], 'LineWidth', 1.5, 'MeshDensity', 2000)
grid on
title('Dial reading vs time', 'Interpreter', 'latex')
xlabel('Time, t (s)', 'Interpreter', 'latex')
ylabel('Indicated mass (kg)', 'Interpreter', 'latex')

% Print Tables
fprintf('k1 each = %.0f N/m, K_B = %.3f N*m/rad, K_T = %.0f N/m\n', k1, K_B, K_T);
fprintf('Lever ratio L_AE/L_AB = %.3f, lever arm gain = %.2f\n', L_AE/L_AB, G_L);
fprintf('J_eq = %.3e kg*m^2, K_eq = %.4f N*m/rad, B_eq = %.3e N*m*s/rad\n', J_eq, K_eq, B_eq);
fprintf('Natural frequency = %.2f Hz\n', sqrt(K_eq/J_eq)/(2*pi));
fprintf('Steady-state theta = %.2f rad (%.0f deg) for %.0f kg\n', ...
        T_in/K_eq, rad2deg(T_in/K_eq), M_person);
fprintf('Mass for one full dial turn = %.0f kg\n', 2*pi*K_eq/(r_E*g));