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

% joint values
% sore the positions for plotting
new_B_x(1) = B(1);
new_B_y(1) = B(2);
new_C_x(1) = C(1);
new_C_y(1) = C(2);
new_E_x(1) = E(1);
new_E_y(1) = E(2);
new_F_x(1) = F(1);
new_F_y(1) = F(2);


% Define the lengths of the bars
lAB = norm(B - A);
lBC = norm(C - B);
lCD = norm(D - C);
lBE = norm(E - B);
lEF = norm(F - E);
lFG = norm(G - F);
lCE = norm (E - C);

% Weight of Each Link
WAB = [0 -1 0];
WBEC = [0 -1 0];
WCD = [0 -1 0];
WEF = [0 -1 0];
WFG = [0 -1 0];

% Center of mass of each link
S1 = (A+B)/2;
S2 = (B+C+E)/3; % Not accurate
S3 = (C+D)/2;
S4 = (E+F)/2;
S5 = (F+G)/2;

syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
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

eqn2 = cross(A-S1, ForceA) + cross(B-S1, ForceB) + InputTorque == 0;

% Equations for Link BEC
% Sum of Forces = 0
% -Fb + Fc + Fe + WBEC = 0
eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;


 % Sum of Moments = 0
% With respect to the center of mass of link BEC

eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == 0;

% -Fc + Fd + WCD = 0
% Equations for Link CD
% Sum of Forces = 0 for Link CD
eqn5 = -ForceC + ForceD + WCD == 0;

% Sum of Moments = 0 for link CD
% Sum of moments = 0 with respect to center of mass of link CD
% S3C x - FC + S3D x FD = 0
% Sum of Moments = 0 for Link CD
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;

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
eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == 0;

% Soving the 10 equations

eqnMatric = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];

StaticSolution = solve (eqnMatric, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);
Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);
Input_Torque = double(StaticSolution.Tin);

disp('Force A:');
disp([Force_Ax, Force_Ay]);
disp('Force B:');
disp([Force_Bx, Force_Bx]);
disp('Force C:');
disp([Force_Cx, Force_Cy]);
disp('Force D:');
disp([Force_Dx, Force_Dy]);
disp('Force E:');
disp([Force_Ex, Force_Ey]);
disp('Force F:');
disp([Force_Fx, Force_Fy]);
disp('Force G:');
disp([Force_Gx, Force_Gy]);
disp('Input Torque:');
disp(Input_Torque);

% Angular Velocity Calculations
% Loop ABCDA

syms wBEC wCD 
omega_AB = [0 0 1];
omega_BEC = [0 0 wBEC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB, B-A) + cross(omega_BEC, C-B) + cross(omega_CD, D-C) == 0;

loopsolution = solve(eqn11, [wBEC wCD]);

% Angular velocity results
angularVelocity_BEC = double(loopsolution.wBEC);
angularVelocity_CD = double(loopsolution.wCD);

omegaBEC = [0 0 angularVelocity_BEC];
omegaCD = [0 0 angularVelocity_CD]; 

syms  wEF wFG 
omega_EF= [0 0 wEF];
omega_FG = [0 0 wFG];

eqn12 = cross(omegaCD, C-D) + cross(omegaBEC, E-C) + cross(omega_EF, F-E) + cross(omega_FG, G-F) == 0

loop2Solution = solve(eqn12, [wEF wFG]);

% Angular velocity results
angularVelocity_EF = double(loop2Solution.wEF)
angularVelocity_FG = double(loop2Solution.wFG)


% Angular Acceleration Calculations
% Define angular accelerations
syms aBEC aCD

alpha_AB = [0 0 0];
alpha_BEC = [0 0 aBEC];
alpha_CD = [0 0 aCD];

a_B_A = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));
a_C_B = cross(alpha_BEC, C-B) + cross(omegaBEC, cross(omegaBEC, C-B));
a_D_C = cross(alpha_CD, D-C) + cross(omegaCD, cross(omegaCD, D-C));

eqn13 = a_B_A + a_C_B + a_D_C == 0;

loop1AccSolution = solve(eqn13, [aBEC aCD]);

% Angular acceleration results
alphaBEC = double(loop1AccSolution.aBEC)
alphaCD = double(loop1AccSolution.aCD)

alphaBEC_vector = [0 0 alphaBEC];
alphaCD_vector = [0 0 alphaCD];

syms aEF aFG

alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

% a_C_D + a_E_C +a_F_E + a_G_F = 0
a_C_D = cross(alphaCD_vector, C-D) + cross(omegaCD, cross(omegaCD, C-D));

a_E_C = cross(alphaBEC_vector, E-C) + cross(omegaBEC, cross(omegaBEC, E-C));

angVel_EF = [0 0 angularVelocity_EF];
angVel_FG = [0 0 angularVelocity_FG];

a_F_E = cross(alpha_EF, F-E) + cross(angVel_EF, cross(angVel_EF, F-E));
a_G_F = cross(alpha_FG, G-F) + cross(angVel_FG, cross(angVel_FG, G-F));

eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

loop2AccSolution = solve(eqn14, [aEF aFG]);

% Angular acceleration results
alphaEF = double(loop2AccSolution.aEF)
alphaFG = double(loop2AccSolution.aFG)

% Velocity at a joint

vB_A = cross(omega_AB, B-A);

% VE_A = V_E_B + V_B_A;

v_E_B = cross(omegaBEC, E-B);

vE_A = v_E_B + vB_A;

%V_S4/G = V_S4_F + V_F_G

V_S4_F = cross(angVel_EF, S4-F);
V_F_G = cross(angVel_FG, F-G);

vS4_G = V_S4_F + V_F_G 

% Velocity of Joint C
%VC_D = omegaCD x C-D

vC_D = cross(omegaCD, C-D);

% Velocity of JOint F
% VF_E = VE_A + omegaEF x F-E
vF_E = vE_A + cross(angVel_EF, F-E);

%Velocity of fixed joints

vA = [0 0 0];
vD = [0 0 0];
vG = [0 0 0];

%Velocity at Center of Mass of Link AB
% VS1_A = omegaAB x S1-A

vS1_A = cross(omega_AB, S1-A);

% Velocity at Center of Mass of Link BEC
% VS2_A = VB_A + VS2_B

vS2_B = cross(omegaBEC, S2-B);

vS2_A = vB_A + vS2_B;

% Velocity at Center of Mass of Link CD

% VS3_D = omegaCD x S3-D

vS3_D = cross(omegaCD, S3-D);

% Velocity at Center of Mass of Link EF

% VS4_G = VS4_F + VF_G

V_S4_F = cross(angVel_EF, S4-F);

V_F_G = cross(angVel_FG, F-G);

vS4_G = V_S4_F + V_F_G;

% Velocity at Center of Mass of Link FG

% VS5_G = omegaFG x S5-G

vS5_G = cross(angVel_FG, S5-G);


% Acceleration at Joint B
% aB_A = alphaAB x B-A + omegaAB x (omegaAB x B-A)

aB_A = cross(alpha_AB, B-A) + cross(omega_AB, cross(omega_AB, B-A));

% Acceleration at Joint C

% aC_D = alphaCD x C-D + omegaCD x (omegaCD x C-D)

aC_D = cross(alphaCD_vector, C-D) + cross(omegaCD, cross(omegaCD, C-D));

% Acceleration at Joint E

% aE_A = aB_A + aE_B
aE_B = cross(alphaBEC_vector, E-B) + cross(omegaBEC, cross(omegaBEC, E-B));

aE_A = aB_A + aE_B;

% Acceleration at Joint F

% aF_E = aE_A + alphaEF x F-E + omegaEF x (omegaEF x F-E)

alphaEF_vector = [0 0 alphaEF];

alphaFG_vector = [0 0 alphaFG];

aF_E = aE_A + cross(alphaEF_vector, F-E) + cross(angVel_EF, cross(angVel_EF, F-E));

% Acceleration at fixed joints

aA = [0 0 0];

aD = [0 0 0];

aG = [0 0 0];

% Acceleration at Center of Mass of Link AB

% aS1_A = alphaAB x S1-A + omegaAB x (omegaAB x S1-A)

aS1_A = cross(alpha_AB, S1-A) + cross(omega_AB, cross(omega_AB, S1-A));

% Acceleration at Center of Mass of Link BEC

% aS2_A = aB_A + aS2_B
aS2_B = cross(alphaBEC_vector, S2-B) + cross(omegaBEC, cross(omegaBEC, S2-B));

aS2_A = aB_A + aS2_B;

% Acceleration at Center of Mass of Link CD

% aS3_D = alphaCD x S3-D + omegaCD x (omegaCD x S3-D)

aS3_D = cross(alphaCD_vector, S3-D) + cross(omegaCD, cross(omegaCD, S3-D));

% Acceleration at Center of Mass of Link EF

% aS4_E = aE_A + alphaEF x S4-E + omegaEF x (omegaEF x S4-E)

aS4_E = aE_A + cross(alphaEF_vector, S4-E) + cross(angVel_EF, cross(angVel_EF, S4-E));

% Acceleration at Center of Mass of Link FG

% aS5_G = alphaFG x S5-G + omegaFG x (omegaFG x S5-G)

aS5_G = cross(alphaFG_vector, S5-G) + cross(angVel_FG, cross(angVel_FG, S5-G));


% Display Angular Velocities

disp('Angular Velocity AB:');

disp(omega_AB(3));
disp('Angular Velocity BEC:');

disp(angularVelocity_BEC);

disp('Angular Velocity CD:');

disp(angularVelocity_CD);

disp('Angular Velocity EF:');

disp(angularVelocity_EF);

disp('Angular Velocity FG:');

disp(angularVelocity_FG);


% Display Angular Accelerations

disp('Angular Acceleration AB:');

disp(alpha_AB(3));

disp('Angular Acceleration BEC:');

disp(alphaBEC);

disp('Angular Acceleration CD:');

disp(alphaCD);

disp('Angular Acceleration EF:');
disp(alphaEF);

disp('Angular Acceleration FG:');

disp(alphaFG);


% Display Joint Velocities

disp('Velocity Joint A:');

disp(vA);

disp('Velocity Joint B:');

disp(vB_A);

disp('Velocity Joint C:');

disp(vC_D);

disp('Velocity Joint D:');
disp(vD);

disp('Velocity Joint E:');

disp(vE_A);

disp('Velocity Joint F:');

disp(vF_E);

disp('Velocity Joint G:');

disp(vG);


% Display Joint Accelerations

disp('Acceleration Joint A:');

disp(aA);

disp('Acceleration Joint B:');

disp(aB_A);

disp('Acceleration Joint C:');

disp(aC_D);

disp('Acceleration Joint D:');

disp(aD);
disp('Acceleration Joint E:');

disp(aE_A);

disp('Acceleration Joint F:');

disp(aF_E);

disp('Acceleration Joint G:');

disp(aG);


% Display Center of Mass Velocities

disp('Velocity Center of Mass Link AB:');

disp(vS1_A);

disp('Velocity Center of Mass Link BEC:');

disp(vS2_A);

disp('Velocity Center of Mass Link CD:');

disp(vS3_D);

disp('Velocity Center of Mass Link EF:');

disp(vS4_G);

disp('Velocity Center of Mass Link FG:');

disp(vS5_G);
% Display Center of Mass Accelerations

disp('Acceleration Center of Mass Link AB:');

disp(aS1_A);

disp('Acceleration Center of Mass Link BEC:');

disp(aS2_A);

disp('Acceleration Center of Mass Link CD:');

disp(aS3_D);

disp('Acceleration Center of Mass Link EF:');

disp(aS4_E);

disp('Acceleration Center of Mass Link FG:');

disp(aS5_G);

% Newton's Second Law Implementation

MassAB = 1;
MassBEC = 1;
MassCD = 1;
MassEF = 1;
MassFG = 1;

% Mass moment of inertia
J_AB = 1;
J_BEC = 1;
J_CD = 1;
J_EF = 1;
J_FG = 1;

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

% Define forces

NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NInputTorque = [0 0 NTin];

% Equations for Link AB
% Sum of forces

eqn15 = NForceA + NForceB + WAB == MassAB * aS1_A;

% Sum of moments = 0
eqn16 = cross(A-S1, NForceA) + cross(B-S1, NForceB) + NInputTorque == J_AB * alpha_AB;

% Equations for Link BEC
% Sum of forces for Link BEC
eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBEC * aS2_A;

% Sum of moments = 0
eqn18 = cross(B-S2, -NForceB) + cross(C-S2, NForceC) + cross(E-S2, NForceE) == J_BEC * alphaBEC_vector;

% Equations for Link CD
% Sum of forces
eqn19 = -NForceC + NForceD + WCD == MassCD * aS3_D;

% SUm of moments
eqn20 = cross(C-S3, -NForceC) + cross(D-S3, NForceD) == J_CD * alphaCD_vector;

% Equations for Link EF
% Sum of forces
eqn21 = -NForceE + NForceF + WEF == MassEF * aS5_G;

% Sum of moments
eqn22 = cross(E-S4, -NForceE) + cross(F-S4, NForceF) == J_EF * [0 0 alphaEF];

% Equations for Link FG
%Sum of forces
eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG * aS5_G;

% Sum of moments 
eqn24 = cross(F-S5, -NForceF) + cross(G-S5, NForceG) == J_FG * [0 0 alphaFG];

%Solving equations
NeqnMatrix = [eqn15 eqn16 eqn17 eqn18 eqn19 eqn20 eqn21 eqn22 eqn23 eqn24];
DynamicSolution = solve(NeqnMatrix, [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin]);

% Extracting the forces from the solution
NForceAx_sol = DynamicSolution.NFAx;
NForceAy_sol = DynamicSolution.NFAy;
NForceBx_sol = DynamicSolution.NFBx;
NForceBy_sol = DynamicSolution.NFBy;
NForceCx_sol = DynamicSolution.NFCx;
NForceCy_sol = DynamicSolution.NFCy;
NForceDx_sol = DynamicSolution.NFDx;
NForceDy_sol = DynamicSolution.NFDy;
NForceEx_sol = DynamicSolution.NFEx;
NForceEy_sol = DynamicSolution.NFEy;
NForceFx_sol = DynamicSolution.NFFx;
NForceFy_sol = DynamicSolution.NFFy;
NForceGx_sol = DynamicSolution.NFGx;
NForceGy_sol = DynamicSolution.NFGy;
NTin_sol = DynamicSolution.NTin;

% Circle Intersection Technique
% Joint coordinates have been defined
% Length of links also defined

% Compute initial angle of input link

initial_theta = atan2(B(2) - A(2), B(1) - A(1));

if (initial_theta < 0)
    inputAngle = 2*pi + initial_theta;
else
    inputAngle = initial_theta;
end


for theta = 1:1:360

    % Find new position of joint B

    B_new = A + [lAB*cos(inputAngle + deg2rad(theta)) lAB*sin(inputAngle + deg2rad(theta)) 0];

    % Find new position of C
    % With B_new as center, BC as radius
    % With D as center and DC as radius

    [Cx, Cy] = circcirc (B_new(1), B_new(2), lBC, D(1), D(2), lCD);

    % Checking if there is a NaN (not a number)

    circIntersect_x_C = any(isnan(vpa(Cx)));
    circIntersect_y_C = any(isnan (vpa(Cy)));
    
    if circIntersect_x_C ==0 && circIntersect_y_C ==0
        C_1 = [Cx(1) Cy(1) 0];
        C_2 = [Cx(2) Cy(2) 0];

        % Distance to determine whether C_1 or C_2 is correct

        dist1 = norm(C_1 - C);
        dist2 = norm(C_2 - C);

        if (dist1<dist2)
            C_new = vpa(C_1);

        else

            C_new = vpa(C_2);
            
        end
   
  % New position of Joint E using B_new and C_new
    [Ex, Ey] = circcirc (B_new(1), B_new(2), lBE, C_new(1), C_new(2), lCE);

    circIntersect_x_E = any(isnan(vpa(Ex)));
    circIntersect_y_E = any(isnan (vpa(Ey)));
    
    if circIntersect_x_E ==0 && circIntersect_y_E ==0
        E_1 = [Ex(1) Ey(1) 0];
        E_2 = [Ex(2) Ey(2) 0];
        % Distance to determine whether E_1 or E_2 is correct

        dist1 = norm(E_1 - E);
        dist2 = norm(E_2 - E);

        if(dist1<dist2)
            E_new = vpa(E_1);
        else
            E_new = vpa(E_2);
        end

    % New position of Joint F using E_new and G
        [Fx, Fy] = circcirc (E_new(1), E_new(2), lEF, G(1), G(2), lFG);
        
        % Checking if there is a NaN
        circIntersect_x_F = any(isnan(vpa(Fx)));
        circIntersect_y_F = any(isnan (vpa(Fy)));
    
        if circIntersect_x_F == 0 && circIntersect_y_F ==0
            F_1 = [Fx(1) Fy(1) 0];
            F_2 = [Fx(2) Fy(2) 0];
    
        % Distance to determine whether E_1 or E_2 is correct

            distF1 = norm(F_1 - F);
            distF2 = norm(F_2 - F);

            if (distF1<distF2)
                F_new = vpa(F_1);
            else
                F_new = vpa(F_2);
            end

% Store values for plotting
            new_B_x(theta+1) = B_new(1);
            new_B_y(theta+1) = B_new(2);
    
            new_C_x(theta+1) = C_new(1);
            new_C_y(theta+1) = C_new(2);
    
            new_E_x(theta+1) = E_new(1);
            new_E_y(theta+1) = E_new(2);
            
            new_F_x(theta+1) = F_new(1);
            new_F_y(theta+1) = F_new(2);
      
       
            B = B_new ;
            C = C_new ;
            E = E_new ;
            F = F_new ;

            % Static equilibrium code
    
            % Center of mass of each link
            S1 = (A+B)/2;
            S2 = (B+C+E)/3; % Not accurate
            S3 = (C+D)/2;
            S4 = (E+F)/2;
            S5 = (F+G)/2;
                
            syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin
            
            ForceA = [FAx FAy 0];
            ForceB = [FBx FBy 0];
            ForceC = [FCx FCy 0];
            ForceD = [FDx FDy 0];
            ForceE = [FEx FEy 0];
            ForceF = [FFx FFy 0];
            ForceG = [FGx FGy 0];
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
            
            eqn2 = cross(A-S1, ForceA) + cross(B-S1, ForceB) + InputTorque == 0;
            
            % Equations for Link BEC
            % Sum of Forces = 0
            % -Fb + Fc + Fe + WBEC = 0
            eqn3 = -ForceB + ForceC + ForceE + WBEC == 0;
            
            
             % Sum of Moments = 0
            % With respect to the center of mass of link BEC
            
            eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) + cross(E-S2, ForceE) == 0;
            
            % -Fc + Fd + WCD = 0
            % Equations for Link CD
            % Sum of Forces = 0 for Link CD
            eqn5 = -ForceC + ForceD + WCD == 0;
            
            % Sum of Moments = 0 for link CD
            % Sum of moments = 0 with respect to center of mass of link CD
            % S3C x - FC + S3D x FD = 0
            % Sum of Moments = 0 for Link CD
            eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;
            
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
            eqn10 = cross(F-S5, -ForceF) + cross(G-S5, ForceG) == 0;
            
            % Soving the 10 equations
            
            eqnMatric = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10];
            
            StaticSolution = solve (eqnMatric, [FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);
            
            Force_Ax(theta + 1) = double(StaticSolution.FAx);
            Force_Ay(theta + 1)  = double(StaticSolution.FAy);
            Force_Bx(theta + 1)  = double(StaticSolution.FBx);
            Force_By(theta + 1)  = double(StaticSolution.FBy);
            Force_Cx(theta + 1)  = double(StaticSolution.FCx);
            Force_Cy(theta + 1)  = double(StaticSolution.FCy);
            Force_Dx(theta + 1)  = double(StaticSolution.FDx);
            Force_Dy(theta + 1)  = double(StaticSolution.FDy);
            Force_Ex(theta + 1)  = double(StaticSolution.FEx);
            Force_Ey(theta + 1)  = double(StaticSolution.FEy);
            Force_Fx(theta + 1)  = double(StaticSolution.FFx);
            Force_Fy(theta + 1)  = double(StaticSolution.FFy);
            Force_Gx(theta + 1)  = double(StaticSolution.FGx);
            Force_Gy(theta + 1)  = double(StaticSolution.FGy);
            Input_Torque = double(StaticSolution.Tin);
              


            % velocity and accln
    % Angular Velocity Calculations
            % Loop ABCDA

            syms wBEC wCD

            omega_AB = [0 0 1];
            omega_BEC = [0 0 wBEC];
            omega_CD = [0 0 wCD];

            eqn11 = cross(omega_AB, B-A) + cross(omega_BEC, C-B) + cross(omega_CD, D-C) == 0;

            loopsolution = solve([eqn11(1) eqn11(2)], [wBEC wCD]);

            angularVelocity_BEC = double(loopsolution.wBEC);
            angularVelocity_CD = double(loopsolution.wCD);

            omegaBEC = [0 0 angularVelocity_BEC];
            omegaCD = [0 0 angularVelocity_CD];


            syms wEF wFG

            omega_EF = [0 0 wEF];
            omega_FG = [0 0 wFG];

            eqn12 = cross(omegaCD, C-D) + cross(omegaBEC, E-C) + ...
                cross(omega_EF, F-E) + cross(omega_FG, G-F) == 0;

            loop2Solution = solve([eqn12(1) eqn12(2)], [wEF wFG]);

            angularVelocity_EF = double(loop2Solution.wEF);
            angularVelocity_FG = double(loop2Solution.wFG);

            angVel_EF = [0 0 angularVelocity_EF];
            angVel_FG = [0 0 angularVelocity_FG];


            % Angular Acceleration Calculations

            syms aBEC aCD

            alpha_AB = [0 0 0];
            alpha_BEC = [0 0 aBEC];
            alpha_CD = [0 0 aCD];

            a_B_A = cross(alpha_AB, B-A) + ...
                cross(omega_AB, cross(omega_AB, B-A));

            a_C_B = cross(alpha_BEC, C-B) + ...
                cross(omegaBEC, cross(omegaBEC, C-B));

            a_D_C = cross(alpha_CD, D-C) + ...
                cross(omegaCD, cross(omegaCD, D-C));

            eqn13 = a_B_A + a_C_B + a_D_C == 0;

            loop1AccSolution = solve([eqn13(1) eqn13(2)], [aBEC aCD]);

            alphaBEC = double(loop1AccSolution.aBEC);
            alphaCD = double(loop1AccSolution.aCD);

            alphaBEC_vector = [0 0 alphaBEC];
            alphaCD_vector = [0 0 alphaCD];


            syms aEF aFG

            alpha_EF = [0 0 aEF];
            alpha_FG = [0 0 aFG];

            a_C_D = cross(alphaCD_vector, C-D) + ...
                cross(omegaCD, cross(omegaCD, C-D));

            a_E_C = cross(alphaBEC_vector, E-C) + ...
                cross(omegaBEC, cross(omegaBEC, E-C));

            a_F_E = cross(alpha_EF, F-E) + ...
                cross(angVel_EF, cross(angVel_EF, F-E));

            a_G_F = cross(alpha_FG, G-F) + ...
                cross(angVel_FG, cross(angVel_FG, G-F));

            eqn14 = a_C_D + a_E_C + a_F_E + a_G_F == 0;

            loop2AccSolution = solve([eqn14(1) eqn14(2)], [aEF aFG]);

            alphaEF = double(loop2AccSolution.aEF);
            alphaFG = double(loop2AccSolution.aFG);

            alphaEF_vector = [0 0 alphaEF];
            alphaFG_vector = [0 0 alphaFG];


            % Velocity at Joints

            vB_A = cross(omega_AB, B-A);

            vC_D = cross(omegaCD, C-D);

            v_E_B = cross(omegaBEC, E-B);

            vE_A = v_E_B + vB_A;

            vF_E = vE_A + cross(angVel_EF, F-E);


            % Velocity at Center of Mass

            vS1_A = cross(omega_AB, S1-A);

            vS2_B = cross(omegaBEC, S2-B);

            vS2_A = vB_A + vS2_B;

            vS3_D = cross(omegaCD, S3-D);

            vS4_E = vE_A + cross(angVel_EF, S4-E);

            vS5_G = cross(angVel_FG, S5-G);


            % Acceleration at Joints

            aB_A = cross(alpha_AB, B-A) + ...
                cross(omega_AB, cross(omega_AB, B-A));

            aC_D = cross(alphaCD_vector, C-D) + ...
                cross(omegaCD, cross(omegaCD, C-D));

            aE_B = cross(alphaBEC_vector, E-B) + ...
                cross(omegaBEC, cross(omegaBEC, E-B));

            aE_A = aB_A + aE_B;

            aF_E = aE_A + cross(alphaEF_vector, F-E) + ...
                cross(angVel_EF, cross(angVel_EF, F-E));


            % Acceleration at Center of Mass

            aS1_A = cross(alpha_AB, S1-A) + ...
                cross(omega_AB, cross(omega_AB, S1-A));

            aS2_B = cross(alphaBEC_vector, S2-B) + ...
                cross(omegaBEC, cross(omegaBEC, S2-B));

            aS2_A = aB_A + aS2_B;

            aS3_D = cross(alphaCD_vector, S3-D) + ...
                cross(omegaCD, cross(omegaCD, S3-D));

            aS4_E = aE_A + cross(alphaEF_vector, S4-E) + ...
                cross(angVel_EF, cross(angVel_EF, S4-E));

            aS5_G = cross(alphaFG_vector, S5-G) + ...
                cross(angVel_FG, cross(angVel_FG, S5-G));


            % Newton's Second Law

            syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

            NForceA = [NFAx NFAy 0];
            NForceB = [NFBx NFBy 0];
            NForceC = [NFCx NFCy 0];
            NForceD = [NFDx NFDy 0];
            NForceE = [NFEx NFEy 0];
            NForceF = [NFFx NFFy 0];
            NForceG = [NFGx NFGy 0];

            NInputTorque = [0 0 NTin];


            % Equations for Link AB

            eqn15 = NForceA + NForceB + WAB == MassAB*aS1_A;

            eqn16 = cross(A-S1, NForceA) + ...
                cross(B-S1, NForceB) + NInputTorque == J_AB*alpha_AB;


            % Equations for Link BEC

            eqn17 = -NForceB + NForceC + NForceE + WBEC == MassBEC*aS2_A;

            eqn18 = cross(B-S2, -NForceB) + ...
                cross(C-S2, NForceC) + ...
                cross(E-S2, NForceE) == J_BEC*alphaBEC_vector;


            % Equations for Link CD

            eqn19 = -NForceC + NForceD + WCD == MassCD*aS3_D;

            eqn20 = cross(C-S3, -NForceC) + ...
                cross(D-S3, NForceD) == J_CD*alphaCD_vector;


            % Equations for Link EF

            eqn21 = -NForceE + NForceF + WEF == MassEF*aS4_E;

            eqn22 = cross(E-S4, -NForceE) + ...
                cross(F-S4, NForceF) == J_EF*alphaEF_vector;


            % Equations for Link FG

            eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG*aS5_G;

            eqn24 = cross(F-S5, -NForceF) + ...
                cross(G-S5, NForceG) == J_FG*alphaFG_vector;


            % Solving equations

            NeqnMatrix = [eqn15(1) eqn15(2) eqn16(3) ...
                eqn17(1) eqn17(2) eqn18(3) ...
                eqn19(1) eqn19(2) eqn20(3) ...
                eqn21(1) eqn21(2) eqn22(3) ...
                eqn23(1) eqn23(2) eqn24(3)];

            DynamicSolution = solve(NeqnMatrix, ...
                [NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy ...
                NFEx NFEy NFFx NFFy NFGx NFGy NTin]);


            % Extracting the forces from the solution

            NForceAx_sol(theta+1) = double(DynamicSolution.NFAx);
            NForceAy_sol(theta+1) = double(DynamicSolution.NFAy);

            NForceBx_sol(theta+1) = double(DynamicSolution.NFBx);
            NForceBy_sol(theta+1) = double(DynamicSolution.NFBy);

            NForceCx_sol(theta+1) = double(DynamicSolution.NFCx);
            NForceCy_sol(theta+1) = double(DynamicSolution.NFCy);

            NForceDx_sol(theta+1) = double(DynamicSolution.NFDx);
            NForceDy_sol(theta+1) = double(DynamicSolution.NFDy);

            NForceEx_sol(theta+1) = double(DynamicSolution.NFEx);
            NForceEy_sol(theta+1) = double(DynamicSolution.NFEy);

            NForceFx_sol(theta+1) = double(DynamicSolution.NFFx);
            NForceFy_sol(theta+1) = double(DynamicSolution.NFFy);

            NForceGx_sol(theta+1) = double(DynamicSolution.NFGx);
            NForceGy_sol(theta+1) = double(DynamicSolution.NFGy);

            NTin_sol(theta+1) = double(DynamicSolution.NTin);

            % Newton's Second Law

   
            else 
                fprintf('New Position of F cannot be determined at angle: %d degree', theta);
            end 

        else
            fprintf('New Position of E cannot be determined at angle: %d degree', theta);
        end
    
    else
        fprintf('New Position of C cannot be determined at angle: %d degree', theta);
    end

end 


