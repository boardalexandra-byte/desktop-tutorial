% Homework 1 ME 4324

% Set up 
clear; close all; clc

% Define joint position vectors
% Units are in meters

A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.2 0];

% Calculate fixed lengths
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D-C);
lDE = norm(E - D);
lCE = norm(E - C); % Needed for rigid triangle CDE
lEF = norm(F - E);
lGF = norm(F - G);

% Calculate the gripper (H) geometry relative to GF
rF_G = F - G;
magG_F = norm(rF_G);
magF_H = 1.843;
magG_H = magG_F + magF_H; 
thetaH_init = atan2(rF_G(1), rF_G(2)); % Using (x,y) convention

% Center of Mass Geometry setup (Link 3 and 5)
mag_rs3_D = 1.2463; 
mag_rs5_G = 2.2253;

% Center of mass 
% Calculate the center of mass for the gripper geometry
H = G + magG_H * [cos(thetaH_init), sin(thetaH_init), 0]
S1 = (A+B)/2
S2 = (B + C)/2
S3 = (C +E)/2
S4 = (E + F)/2
S5 = (G + F)/2

% Parameters of aisi steel
M_AB = 23.34; J_AB = 0.657;
M_BC = 56.43; J_BC = 9.4718;
M_DE = 97.35; J_DE = 50.2422;
M_EF = 47.00; J_EF = 5.4503;
M_FG = 174.63; J_FG = 285.3794;

% Weight vectors
g = -9.81;
W1 = [0 (M_AB*g) 0];
W2 = [0 (M_BC*g) 0];
W3 = [0 (M_DE*g) 0];
W4 = [0 (M_EF*g) 0];
W5 = [0 (M_FG*g) 0];
W_H = [0 -200 0];

% Initial Input Angle
initial_angle_AB = atan2(B(2)-A(2), B(1)-A(1));

% Initialize Arrays
steps = 361; % 0 to 360 includes 361 steps

% --- Pre-allocate Storage Arrays (For Plotting Requirement) ---
% Angular Velocities (rad/s)
Store_Omega_BC = zeros(1, steps); Store_Omega_CDE = zeros(1, steps);
Store_Omega_EF = zeros(1, steps); Store_Omega_FGH = zeros(1, steps);

% Angular Accelerations (rad/s^2)
Store_Alpha_BC = zeros(1, steps); Store_Alpha_CDE = zeros(1, steps);
Store_Alpha_EF = zeros(1, steps); Store_Alpha_FGH = zeros(1, steps);

% Linear Accelerations of Joints (m/s^2) - Magnitude only for simple plotting
Store_Acc_B = zeros(1, steps); Store_Acc_C = zeros(1, steps);
Store_Acc_E = zeros(1, steps); Store_Acc_F = zeros(1, steps);

% Dynamic Forces (Newtons) - Magnitude only for simple plotting
Store_Force_A = zeros(1, steps); Store_Force_B = zeros(1, steps);
Store_Force_C = zeros(1, steps); Store_Force_D = zeros(1, steps);
Store_Vel_B = zeros(1, steps); Store_Vel_C = zeros(1, steps);
Store_Vel_E = zeros(1, steps); Store_Vel_F = zeros(1, steps); 
   
% Position Traces
newB_x = zeros(1, steps); newB_y = zeros(1, steps);
newC_x = zeros(1, steps); newC_y = zeros(1, steps);
newE_x = zeros(1, steps); newE_y = zeros(1, steps);
newF_x = zeros(1, steps); newF_y = zeros(1, steps);
newH_x = zeros(1, steps); newH_y = zeros(1, steps);

% Torques
Static_Torque = zeros(1, steps);
Newton_Torque = zeros(1, steps);

% Simulation Loop
for theta = 0:1:360

    idx = theta + 1;

    % Position Analysis

    % New B
    angle_curr = initial_angle_AB + deg2rad(theta);
    B_new = A + [lAB*cos(angle_curr) lAB*sin(angle_curr) 0];

    % New C (@ the intersection B-C and D-C)
    [Cx_sol, Cy_sol] = circcirc (B_new(1), B_new(2), lBC, D(1), D(2), lCD);
    
    if any(isnan(Cx_sol)); fprintf('Singularity at theta=%d\n', theta); break; end 
    
    C1 = [Cx_sol(1) Cy_sol(1) 0];
    C2 = [Cx_sol(2) Cy_sol(2) 0];
    
    % Continuity check (closest to previous C)
    if theta == 0
        prev_C = C;
    else
        prev_C = [newC_x(idx-1) newC_y(idx-1) 0];
    end
    
    if norm(C1 - prev_C) < norm(C2 - prev_C)
        C_new = C1;
    else
        C_new = C2;
    end
    
    % New E (Intersection D-E and C-E) - Rigid Triangle
    [Ex_sol, Ey_sol] = circcirc(D(1), D(2), lDE, C_new(1), C_new(2),lCE);

    if any (isnan(Ex_sol)); fprintf('Singularity at theta=%d\n', theta); break;end

    E1 = [Ex_sol(1) Ey_sol(1) 0];
    E2 = [Ex_sol(2) Ey_sol(2) 0];

    if theta == 0
        prev_E = E;
    else
        prev_E = [newE_x(idx - 1) newE_y(idx - 1) 0];
    end

    if norm(E1 - prev_E) < norm(E2 - prev_E)
        E_new = E1;
    else
        E_new = E2;
    end
    
    % New F ( @ intersection EF & GF)
    [Fx_sol, Fy_sol] = circcirc (E_new(1), E_new(2), lEF, G(1), G(2), lGF);
    if any (isnan(Fx_sol)); fprintf('Singularity at theta=%d\n', theta); break;end

    F1 = [Fx_sol(1) Fy_sol(1) 0];
    F2 = [Fx_sol(2) Fy_sol(2) 0];
    
    if theta == 0
        prev_F = F;
    else
        prev_F = [newF_x(idx - 1) newF_y(idx - 1) 0];
    end

    if norm(F1 - prev_F) < norm(F2 - prev_F)
        F_new = F1;
    else
        F_new = F2;
    end

    % New H of the gripper (extension from FG)
    vec_GF = F_new - G;
    unit_GF = vec_GF / norm(vec_GF);
    H_new = G + unit_GF * magG_H;

    % Store positions for plotting
    newB_x(idx) = B_new(1);
    newB_y(idx) = B_new(2);
    newC_x(idx) = C_new(1);
    newC_y(idx) = C_new(2);
    newE_x(idx) = E_new(1);
    newE_y(idx) = E_new(2);
    newF_x(idx) = F_new(1);
    newF_y(idx) = F_new(2);
    H_positions(idx, :) = H_new; % Store H positions

    % Update the current positions to the updated
    B_curr = B_new;
    C_curr = C_new;
    E_curr = E_new;
    F_curr = F_new;
    H_curr = H_new;

    % Update new center of masses
    S1 = (A + B_curr)/2;
    S2 = (B_curr + C_curr)/2;
    
    % Updating S3
    E_D = E_curr - D;
    theta_s3 = atan2(E_D(1), E_D(2));
    S3 = D + [mag_rs3_D*sin(theta_s3) mag_rs3_D*cos(theta_s3) 0];

    S4 = (E_curr + F_curr)/2;
    
    % Updating S5 
    F_G = F_curr - G;
    theta_s5 = atan2(F_G(1), F_G(2));
    S5 = G + [mag_rs5_G * sin(theta_s5) mag_rs5_G * cos(theta_s5) 0];

    % Static Equilibrium Analysis

    syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tinz
   
    % Force vectors
    FA = [FAx FAy 0];
    FB = [FBx FBy 0];
    FC = [FCx FCy 0];
    FD = [FDx FDy 0];
    FE = [FEx FEy 0];
    FF = [FFx FFy 0];
    FG = [FGx FGy 0];

    Tin = [0 0 Tinz];

    % Equations
    eqn1 = FA + FB + W1 ==0;
    eqn2 = cross(A-S1, FA) + cross(B_curr-S1, FB) + Tin == 0;
    eqn3 = FC - FB + W2 == 0;
    eqn4 = cross(C_curr-S2, FC) + cross(B_curr-S2, -FB) == 0;
    eqn5 = FE - FC + FD + W3 ==0;
    eqn6 = cross(E_curr-S3, FE) + cross(C_curr-S3, -FC) + cross(D-S3, FD) == 0;
    eqn7 = -FE + FF + W4 == 0;
    eqn8 = cross(E_curr-S4, -FE) + cross(F_curr-S4, FF) == 0;
    eqn9 = -FF + FG + W5 + W_H == 0;
    eqn10 = cross(F_curr-S5, -FF) + cross(G-S5, FG) + cross(H_curr-S5, W_H) == 0;

    % Solve equations
    eqns = [eqn1; eqn2; eqn3; eqn4; eqn5; eqn6; eqn7; eqn8; eqn9; eqn10];
    vars = [FAx, FAy, FBx, FBy, FCx, FCy, FDx, FDy, FEx, FEy, FFx, FFy, FGx, FGy, Tinz];
    
    sol_static = solve(eqns, vars);
    Static_Torque(idx) = double(sol_static.Tinz);

    % Velocity Analysis
    syms w_bc w_cde w_ef w_fgh
    w_in = [0 0 2.424];

    % Loop 1 (ABCDEFGA)
    eqn11 = cross(w_in, B_curr-A) + cross([0 0 w_bc], C_curr-B_curr) + cross([0 0 w_cde], D-C_curr) == 0;
    sol_v1 = solve(eqn11, [w_bc, w_cde]);
    W_BC = double(sol_v1.w_bc); W_CDE = double(sol_v1.w_cde);

    % Loop 2 (DEFHGD)
    eqn12 = cross([0 0 W_CDE], E_curr-D) + cross([0 0 w_ef], F_curr-E_curr) + cross([0 0 w_fgh], G-F_curr) == 0;
    sol_v2 = solve(eqn12, [w_ef, w_fgh]);
    W_EF = double(sol_v2.w_ef); W_FGH = double(sol_v2.w_fgh);

    % Store velocities
    Store_Omega_BC(idx) = W_BC; Store_Omega_CDE(idx) = W_CDE;
    Store_Omega_EF(idx) = W_EF; Store_Omega_FGH(idx) = W_FGH;
    v_B = cross(w_in, B_curr-A);
    v_C = v_B + cross([0 0 W_BC], C_curr-B_curr);
    v_E = cross([0 0 W_CDE], E_curr-D);
    v_F = v_E + cross([0 0 W_EF], F_curr-E_curr);
    
    Store_Vel_B(idx) = norm(v_B);
    Store_Vel_C(idx) = norm(v_C);
    Store_Vel_E(idx) = norm(v_E);
    Store_Vel_F(idx) = norm(v_F);

    % Acceleration Analysis
    syms a_bc a_cde a_ef a_fgh
    a_in = [0 0 0];

    %Loop 1 
    acc_B = cross(a_in, B_curr-A) + cross(w_in, cross(w_in, B_curr-A));
    acc_CB = cross([0 0 a_bc], C_curr-B_curr) + cross([0 0 W_BC], cross([0 0 W_BC], C_curr-B_curr));
    acc_DC = cross([0 0 a_cde], D-C_curr) + cross([0 0 W_CDE], cross([0 0 W_CDE], D-C_curr));

    eqn13 = acc_B + acc_CB + acc_DC == 0;
    sol_a1 = solve(eqn13, [a_bc, a_cde]);
    A_BC = double(sol_a1.a_bc); A_CDE = double(sol_a1.a_cde);
    
    % Loop 2
    acc_E = cross([0 0 A_CDE], E_curr-D) + cross([0 0 W_CDE], cross([0 0 W_CDE], E_curr-D));
    acc_FE = cross([0 0 a_ef], F_curr-E_curr) + cross([0 0 W_EF], cross([0 0 W_EF], F_curr-E_curr));
    acc_GF = cross([0 0 a_fgh], G-F_curr) + cross([0 0 W_FGH], cross([0 0 W_FGH], G-F_curr));
    
    eqn14 = acc_E + acc_FE + acc_GF == 0;
    sol_a2 = solve(eqn14, [a_ef, a_fgh]);
    A_EF = double(sol_a2.a_ef); A_FGH = double(sol_a2.a_fgh);

    % Store Angular Accelerations
    Store_Alpha_BC(idx) = A_BC; Store_Alpha_CDE(idx) = A_CDE;
    Store_Alpha_EF(idx) = A_EF; Store_Alpha_FGH(idx) = A_FGH;

     % Store Linear Accelerations (Joints)
    Store_Acc_B(idx) = norm(acc_B);
    acc_C_val = acc_B + cross([0 0 A_BC], C_curr-B_curr) + cross([0 0 W_BC], cross([0 0 W_BC], C_curr-B_curr));
    Store_Acc_C(idx) = norm(acc_C_val);
    Store_Acc_E(idx) = norm(acc_E);
    
    % Accel F
    acc_F_val = acc_E + cross([0 0 A_EF], F_curr-E_curr) + cross([0 0 W_EF], cross([0 0 W_EF], F_curr-E_curr));
    Store_Acc_F(idx) = norm(acc_F_val);

     % COM Accelerations
    Acc_S1 = cross(a_in, S1-A) + cross(w_in, cross(w_in, S1-A));
    Acc_S2 = acc_B + cross([0 0 A_BC], S2-B_curr) + cross([0 0 W_BC], cross([0 0 W_BC], S2-B_curr));
    Acc_S3 = cross([0 0 A_CDE], S3-D) + cross([0 0 W_CDE], cross([0 0 W_CDE], S3-D));
    Acc_S4 = acc_E + cross([0 0 A_EF], S4-E_curr) + cross([0 0 W_EF], cross([0 0 W_EF], S4-E_curr));
    Acc_S5 = cross([0 0 A_FGH], S5-G) + cross([0 0 W_FGH], cross([0 0 W_FGH], S5-G));
    
    % --- Newton's Second Law ---
    syms N_FAx N_FAy N_FBx N_FBy N_FCx N_FCy N_FDx N_FDy N_FEx N_FEy N_FFx N_FFy N_FGx N_FGy N_Tin
    
    % Dynamic Force Vectors
    N_FA=[N_FAx N_FAy 0]; N_FB=[N_FBx N_FBy 0]; N_FC=[N_FCx N_FCy 0];
    N_FD=[N_FDx N_FDy 0]; N_FE=[N_FEx N_FEy 0]; N_FF=[N_FFx N_FFy 0]; N_FG=[N_FGx N_FGy 0];
    N_Torque = [0 0 N_Tin];
    
    eqn15 = N_FA + N_FB + W1 == M_AB * Acc_S1;
    eqn16 = cross(A-S1, N_FA) + cross(B_curr-S1, N_FB) + N_Torque == J_AB * a_in;
    eqn17 = -N_FB + N_FC + W2 == M_BC * Acc_S2;
    eqn18 = cross(B_curr-S2, -N_FB) + cross(C_curr-S2, N_FC) == J_BC * [0 0 A_BC];
    eqn19 = -N_FC + N_FD + N_FE + W3 == M_DE * Acc_S3;
    eqn20 = cross(C_curr-S3, -N_FC) + cross(D-S3, N_FD) + cross(E_curr-S3, N_FE) == J_DE * [0 0 A_CDE];
    eqn21 = -N_FE + N_FF + W4 == M_EF * Acc_S4;
    eqn22 = cross(E_curr-S4, -N_FE) + cross(F_curr-S4, N_FF) == J_EF * [0 0 A_EF];
    eqn23 = -N_FF + N_FG + W5 + W_H == M_FG * Acc_S5;
    eqn24 = cross(F_curr-S5, -N_FF) + cross(G-S5, N_FG) + cross(H_curr-S5, W_H) == J_FG * [0 0 A_FGH];
    
    N_eqns = [eqn15; eqn16; eqn17; eqn18; eqn19; eqn20; eqn21; eqn22; eqn23; eqn24];
    N_vars = [N_FAx, N_FAy, N_FBx, N_FBy, N_FCx, N_FCy, N_FDx, N_FDy, N_FEx, N_FEy, N_FFx, N_FFy, N_FGx, N_FGy, N_Tin];
    
    solution_Newton = solve(N_eqns, N_vars);
    
    % Store Force Magnitudes
    FA_vec = [double(solution_Newton.N_FAx) double(solution_Newton.N_FAy) 0];
    FB_vec = [double(solution_Newton.N_FBx) double(solution_Newton.N_FBy) 0];
    FC_vec = [double(solution_Newton.N_FCx) double(solution_Newton.N_FCy) 0];
    FD_vec = [double(solution_Newton.N_FDx) double(solution_Newton.N_FDy) 0];
    FE_vec = [double(solution_Newton.N_FEx) double(solution_Newton.N_FEy) 0];
    FF_vec = [double(solution_Newton.N_FFx) double(solution_Newton.N_FFy) 0];
    FG_vec = [double(solution_Newton.N_FGx) double(solution_Newton.N_FGy) 0];
    
    Store_Force_A(idx) = norm(FA_vec);
    Store_Force_B(idx) = norm(FB_vec);
    Store_Force_C(idx) = norm(FC_vec);
    Store_Force_D(idx) = norm(FD_vec);
    Store_Force_E(idx) = norm(FE_vec);
    Store_Force_F(idx) = norm(FF_vec);
    Store_Force_G(idx) = norm(FG_vec);
    
    Newton_Torque(idx) = double(solution_Newton.N_Tin);


    if theta == 0
        fprintf('\n=======================================\n');
        fprintf('   FIRST POSITION DATA (THETA = 0)   \n');
        fprintf('=======================================\n');
        
        % 1. STATIC EQUILIBRIUM
        fprintf('\n--- 1. STATIC EQUILIBRIUM ---\n');
        fprintf('Static Input Torque: %.4f Nm\n', double(sol_static.Tinz));
        fprintf('Static Force A: %.4f N\n', norm([double(sol_static.FAx) double(sol_static.FAy) 0]));
        fprintf('Static Force B: %.4f N\n', norm([double(sol_static.FBx) double(sol_static.FBy) 0]));
        fprintf('Static Force C: %.4f N\n', norm([double(sol_static.FCx) double(sol_static.FCy) 0]));
        fprintf('Static Force D: %.4f N\n', norm([double(sol_static.FDx) double(sol_static.FDy) 0]));
        fprintf('Static Force E: %.4f N\n', norm([double(sol_static.FEx) double(sol_static.FEy) 0]));
        fprintf('Static Force F: %.4f N\n', norm([double(sol_static.FFx) double(sol_static.FFy) 0]));
        fprintf('Static Force G: %.4f N\n', norm([double(sol_static.FGx) double(sol_static.FGy) 0]));
        
        % 2. NEWTON-EULER (DYNAMIC)
        fprintf('\n--- 2. NEWTON-EULER (DYNAMIC) ---\n');
        fprintf('Dynamic Input Torque: %.4f Nm\n', Newton_Torque(idx));
        fprintf('Dynamic Force A: %.4f N\n', norm(FA_vec));
        fprintf('Dynamic Force B: %.4f N\n', norm(FB_vec));
        fprintf('Dynamic Force C: %.4f N\n', norm(FC_vec));
        fprintf('Dynamic Force D: %.4f N\n', norm(FD_vec));
        fprintf('Dynamic Force E: %.4f N\n', norm(FE_vec));
        fprintf('Dynamic Force F: %.4f N\n', norm(FF_vec));
        fprintf('Dynamic Force G: %.4f N\n', norm(FG_vec));

        % 3. KINEMATICS - ANGULAR
        fprintf('\n--- 3. ANGULAR KINEMATICS ---\n');
        fprintf('Omega (rad/s):   BC=%.4f, CDE=%.4f, EF=%.4f, FGH=%.4f\n', W_BC, W_CDE, W_EF, W_FGH);
        fprintf('Alpha (rad/s^2): BC=%.4f, CDE=%.4f, EF=%.4f, FGH=%.4f\n', A_BC, A_CDE, A_EF, A_FGH);
        
        % 4. KINEMATICS - JOINT LINEAR VELOCITIES
        % Calculate linear velocities just for reporting
        v_B = cross(w_in, B_curr-A);
        v_C = v_B + cross([0 0 W_BC], C_curr-B_curr);
        v_E = cross([0 0 W_CDE], E_curr-D);
        v_F = v_E + cross([0 0 W_EF], F_curr-E_curr);
        
        fprintf('\n--- 4. JOINT LINEAR VELOCITIES (m/s) ---\n');
        fprintf('Mag Vel B: %.4f\n', norm(v_B));
        fprintf('Mag Vel C: %.4f\n', norm(v_C));
        fprintf('Mag Vel E: %.4f\n', norm(v_E));
        fprintf('Mag Vel F: %.4f\n', norm(v_F));

        % 5. KINEMATICS - JOINT LINEAR ACCELERATIONS
        fprintf('\n--- 5. JOINT LINEAR ACCELERATIONS (m/s^2) ---\n');
        fprintf('Mag Acc B: %.4f\n', norm(acc_B));
        fprintf('Mag Acc C: %.4f\n', norm(acc_C_val));
        fprintf('Mag Acc E: %.4f\n', norm(acc_E));
        fprintf('Mag Acc F: %.4f\n', norm(acc_F_val));

        % 6. KINEMATICS - MASS CENTER ACCELERATIONS
        fprintf('\n--- 6. MASS CENTER ACCELERATIONS (m/s^2) ---\n');
        fprintf('Mag Acc S1 (Link AB):  %.4f\n', norm(Acc_S1));
        fprintf('Mag Acc S2 (Link BC):  %.4f\n', norm(Acc_S2));
        fprintf('Mag Acc S3 (Link DE): %.4f\n', norm(Acc_S3));
        fprintf('Mag Acc S4 (Link EF):  %.4f\n', norm(Acc_S4));
        fprintf('Mag Acc S5 (Link GH): %.4f\n', norm(Acc_S5));
        
        fprintf('=======================================\n\n');
    end
    
    % Optional: Progress indicator
    if mod(theta, 15) == 0
        fprintf('Calculated Angle: %d / 360\n', theta);
    end

end %for loop

% 3. Plotting Results (Full 2-Figure Set)

angles_vec = 0:1:360; 

% --- FIGURE 1: Trajectories and Torques ---
figure('Name', 'Mechanism Analysis - Motion & Torque', 'Color', 'w');

subplot(3,2,1); plot(newB_x, newB_y); title('Path of B'); axis equal; grid on; xlabel('X(m)'); ylabel('Y(m)');
subplot(3,2,2); plot(newC_x, newC_y); title('Path of C'); axis equal; grid on; xlabel('X(m)'); ylabel('Y(m)');
subplot(3,2,3); plot(newE_x, newE_y); title('Path of E'); axis equal; grid on; xlabel('X(m)'); ylabel('Y(m)');
subplot(3,2,4); plot(newF_x, newF_y); title('Path of F'); axis equal; grid on; xlabel('X(m)'); ylabel('Y(m)');

subplot(3,2,5); plot(angles_vec, Static_Torque); 
title('Static Torque'); grid on; xlabel('Angle (deg)'); ylabel('Torque (Nm)');

subplot(3,2,6); plot(angles_vec, Newton_Torque); 
title('Dynamic Torque'); grid on; xlabel('Angle (deg)'); ylabel('Torque (Nm)');

% --- FIGURE 2: Kinematics & Kinetics (Velocities, Accels, Forces) ---
figure('Name', 'Mechanism Analysis - Kinematics & Kinetics', 'Color', 'w');

% Angular Velocities
subplot(3,2,1); 
hold on; 
grid on;

% Input link AB
plot(angles_vec, Store_Omega_BC, 'b'); 
plot(angles_vec, Store_Omega_CDE, 'r');
plot(angles_vec, Store_Omega_EF, 'm');
plot(angles_vec, Store_Omega_FGH, 'g');

legend('Link BC', 'Link CDE', 'Link EF', 'Link FGH');
title('Angular Velocities'); ylabel('rad/s');

% Angular Accelerations
subplot(3,2,2); hold on; grid on;
plot(angles_vec, Store_Alpha_BC, 'b'); 
plot(angles_vec, Store_Alpha_CDE, 'r');
plot(angles_vec, Store_Alpha_EF, 'm');
plot(angles_vec, Store_Alpha_FGH, 'g');

legend('Link BC', 'Link CDE', 'Link EF', 'Link FGH'); 

title('Angular Accelerations'); ylabel('rad/s^2');

% Linear Accelerations (Joints)
subplot (3, 2, 3);
hold on;
grid on;
plot(angles_vec, Store_Acc_B, 'b'); 
plot(angles_vec, Store_Acc_C, 'g');
plot(angles_vec, Store_Acc_E, 'r');
plot(angles_vec, Store_Acc_F, 'm');

legend('Joint B', 'Joint C', 'Joint E', 'Joint F');
title('Joint Linear Accel (Mag)'); ylabel('m/s^2');

% Dynamic Forces (Reactions)
subplot(3,2,4); hold on; grid on;
plot(angles_vec, Store_Force_A, 'k'); 
plot(angles_vec, Store_Force_D, 'm');
plot(angles_vec, Store_Force_G, 'b');

legend('Reaction @ A', 'Reaction @ D', 'Reaction @ G');
title('Ground Reaction Forces'); ylabel('Force (N)');

% Dynamic Forces (Joints)
subplot(3,2,5); hold on; grid on;
plot(angles_vec, Store_Force_B, 'b'); 
plot(angles_vec, Store_Force_C, 'g');
plot(angles_vec, Store_Force_E, 'r');
plot(angles_vec, Store_Force_F, 'm');

legend('Joint B', 'Joint C', 'Joint E', 'Joint F');
title('Internal Joint Forces'); ylabel('Force (N)');
% Linear Velocities (Joints)
subplot(3,2,6); hold on; grid on;
plot(angles_vec, Store_Vel_B, 'k--');
plot(angles_vec, Store_Vel_C, 'b');
plot(angles_vec, Store_Vel_E, 'r');
plot(angles_vec, Store_Vel_F, 'g');

legend('Joint B', 'Joint C', 'Joint E', 'Joint F');
title('Joint Linear Velocities (Mag)');
ylabel('m/s');
xlabel('Angle (deg)');
disp('Analysis Complete. All plots generated.');


% 3. Plotting Results 

figure('Name', 'Mechanism Analysis Results', 'NumberTitle', 'off');

% Plot 1: Joint Traces (Mechanism Trajectory)
subplot(2,1,1);
hold on; axis equal; grid on;

% Plot the paths using the arrays generated in the loop
plot(newB_x, newB_y, 'b', 'LineWidth', 1.5);
plot(newC_x, newC_y, 'g', 'LineWidth', 1.5);
plot(newE_x, newE_y, 'm', 'LineWidth', 1.5);
plot(newF_x, newF_y, 'k', 'LineWidth', 1.5);
% Calculate Gripper (H) trace if not explicitly stored in loop, 
% or if you added newH_x/newH_y to the loop, use those. 
% (Assuming H was calculated but maybe not stored in arrays in the slow loop example,
%  if you added newH_x to the loop, use: plot(newH_x, newH_y, ...))

% Add Ground points for context
scatter([A(1) D(1) G(1)], [A(2) D(2) G(2)], 100, 'k', 'filled'); 

legend('Path B', 'Path C', 'Path E', 'Path F', 'Ground Joints');
title('Joint Trajectories over 360 Deg Rotation');
xlabel('X Position (m)'); 
ylabel('Y Position (m)');

% --- Plot 2: Dynamic Torque ---
subplot(2,1,2);
hold on; grid on;

% Create an x-axis vector for degrees
angles = 1:length(Newton_Torque);

plot(angles, Newton_Torque, 'r', 'LineWidth', 1.5);
yline(0, 'k--'); % Zero reference line

title('Required Input Torque (Dynamic) - Newton-Euler Method');
xlabel('Input Angle (deg)'); 
ylabel('Torque (Nm)');
xlim([0 360]);

disp('Analysis Complete.');