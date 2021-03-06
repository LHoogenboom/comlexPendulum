%%% INPUTS
[v, l_1, l_2, Mc, m, M1, M2, J1, J2, g, ks, b1, Fbr, bj] = inputs;
%%% Generalized coordinates
% - x: position of the cart 
% - theta_1: rotational angle of rod 1 with respect to the vertical
% - theta_2: rotational angle of rod 2, with respect ot the vertical
syms x(t) theta_1(t)  theta_2(t)
x = str2sym('x(t)');
theta_1 = str2sym('theta_1(t)');
theta_2 = str2sym('theta_2(t)');
%% Defining the coordinates of the masses
%%% m1     
x_1 = x + l_1*sin(theta_1);
y_1 = -l_1*cos(theta_1);
 
%%% m2 
x_2 = x + l_1*sin(theta_2);
y_2 = -l_2*cos(theta_2);

%%% m3
d_1 = l_1*cos((theta_2 - theta_1)/2);
d_2 = sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2);
L = d_1 + d_2;
theta_3 = (theta_1 + theta_2)/2;
x_3 = x + (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2))*sin((theta_1 + theta_2)/2); 
y_3 = -(l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2))*cos((theta_1 + theta_2)/2);

%%% Rods' centers of gravity 
%%% rod 1
x_cg1 = x + (l_1/2)*sin(theta_1);
y_cg1 = -(l_1/2)*cos(theta_1);

%%% rod 2
x_cg2 = x + (l_1/2)*sin(theta_2);
y_cg2 = - (l_1/2)*cos(theta_2);

%%% rod 3 
AB = sqrt((l_2/2)^2 - (sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2 + (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2); 
nm = l_1^2 + (l_2/2)^2 - (sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2 + (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2 - (l_2/2)^2;
dn = 2*l_1*sqrt((l_2/2)^2 - (sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2 + (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2);
beta_2 = acos((l_1^2 + (l_2/2)^2 - (sqrt(l_2^2 - l_1^2*(sin((thetfetaa_2 - theta_1)/2))^2)/2)^2 + (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2 - (l_2/2)^2)/(2*l_1*sqrt((l_2/2)^2 - (sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2 + (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2)/2)^2)));
theta_5 = theta_1 + beta_2;
x_cg3 = x + AB*sin(theta_5);
y_cg3 = -AB*cos(theta_5);

%%% rod 4 
AC = AB;
theta_6 = theta_2 - beta_2;
x_cg4 = x + AC*sin(theta_6); 
y_cg4 = -AC*cos(theta_6); 

%% Derive the linear velocities by finding the time derivatives
%%% Mc
x_dot = diff(x,t);
pretty(x_dot)

%%% m1 
x_dot_1 = diff(x_1,t);
pretty(x_dot_1) % show in comand window
y_dot_1 = diff(y_1,t); 
pretty(y_dot_1)
v1_sq = x_dot_1^2 + y_dot_1^2; %squared velocity of m1

%%% m2: connected with a spring to m1
x_dot_2 = diff(x_2,t);
pretty(x_dot_2) % show in comand window
y_dot_2 = diff(y_2,t); 
pretty(y_dot_2)
v2_sq = x_dot_2^2 + y_dot_2^2; %squared velocity of m2

%%% m3 
x_dot_3 = diff(x_3,t);
pretty(x_dot_3) % show in comand window
y_dot_3 = diff(y_3,t); 
pretty(y_dot_3)
v3_sq = x_dot_3^2 + y_dot_3^2; %squared velocity of m3

%%% center of gravity(cg1) 
x_dot_cg1 = diff(x_cg1,t);
pretty(x_dot_cg1)
y_dot_cg1 = diff(y_cg1,t);
pretty(y_dot_cg1)
vcg1_sq = x_dot_cg1^2 + y_dot_cg1^2;

%%% cg2 
x_dot_cg2 = diff(x_cg2,t);
pretty(x_dot_cg2)
y_dot_cg2 = diff(y_cg2,t);
pretty(y_dot_cg2)
vcg2_sq = x_dot_cg2^2 + y_dot_cg2^2;

%%% cg3 
x_dot_cg3 = diff(x_cg3,t);
pretty(x_dot_cg3)
y_dot_cg3 = diff(y_cg3,t);
pretty(y_dot_cg3)
vcg3_sq = x_dot_cg3^2 + y_dot_cg3^2;

%%% cg4 
x_dot_cg4 = diff(x_cg4,t);
pretty(x_dot_cg4)
y_dot_cg4 = diff(y_cg4,t);
pretty(y_dot_cg4)
vcg4_sq = x_dot_cg4^2 + y_dot_cg4^2;
%% Angular velocities of the rods
%%% rod 1
theta_1_dot = diff(theta_1,t); 
pretty(theta_1_dot)

%%% rod 2 
theta_2_dot = diff(theta_2,t);
pretty(theta_2_dot)

%%% rod 3 
eta_1 = pi - acos((l_2^2 + l_1^2 - (l_1*cos((theta_2 - theta_1)/2) + sqrt(l_2^2 - l_1^2*(sin((theta_2 - theta_1)/2))^2))^2)/(2*l_2*l_1));
eta_1_dot = diff(eta_1,t); 
pretty(eta_1)

%%% rod 4 
eta_2 = eta_1; 
eta_2_dot = diff(eta_2,t);
pretty(eta_2); 

%% Kinetic Energy of the moving parts 
%%% Transportational(?) movement
%%% K_Mc 
K_Mc = (1/2)*Mc*x_dot^2; 

%%% K_m1
K_m1 = (1/2)*m*v1_sq;

%%% K_m2 
K_m2 = (1/2)*m*v2_sq;

%%% K_m3 
K_m3 = (1/2)*m*v3_sq;

%%% K_cg1
K_cg1 = (1/2)*M1*vcg1_sq;

%%% K_cg2
K_cg2 = (1/2)*M1*vcg2_sq;

%%% K_cg3
K_cg3 = (1/2)*M2*vcg3_sq;

%%% K_cg4
K_cg4 = (1/2)*M2*vcg4_sq;

%%% Rotational movement 
%%% rod 1 
K_rod1 = (1/2)*J1*theta_1_dot^2;

%%% rod 2 
K_rod2 = (1/2)*J1*theta_2_dot^2;

%%% rod 3 
K_rod3 = (1/2)*J2*eta_1_dot^2;

%%% rod 4 
K_rod4 = (1/2)*J2*eta_2_dot^2;

%%% Total Kinetic energy 
K_total = K_Mc + K_m1 + K_m2 + K_m3 + K_cg1 + K_cg2 + K_cg3 + K_cg4 + K_rod1 + K_rod2 + K_rod3 + K_rod4; 

%% Potential Energy 
%%% Potential energy due to height difference
%%% P_Mc 
P_Mc = 0;

%%% P_m1
P_m1 = m*g*y_1;

%%% P_m2 
P_m2 = m*g*y_2;

%%% P_m3 
P_m3 = m*g*y_3;

%%% P_cg1
P_cg1 = M1*g*y_cg1;

%%% P_cg2
P_cg2 = M1*g*y_cg2;

%%% P_cg3
P_cg3 = M2*g*y_cg3;

%%% P_cg4
P_cg4 = M2*g*y_cg4;

%%%Potential energy due to the spring 
dx = x_2 - x_1; % spring deformation, assume deformation only over x-direction
P_spring = (1/2)*ks*dx^2 ;

%%% Total Potential Energy 
P_total = P_Mc +P_m1 + P_m2 + P_m3 + P_cg1 + P_cg2 + P_cg3 + P_cg4 + P_spring; 


%% Dissipation Energy 
%%% From friction of the ground
V1 = (1/2)*b1*x_dot^2;

%%% Friction at joints of mass1, mass2, mass3
V2 = (1/2)*bj*(x_dot_1^2 + y_dot_1^2 + x_dot_2^2 + y_dot_2^2 + x_dot_3^2 + y_dot_3^2);

%%% Total Dissipation energy 
V_total = V1 + V2; 

%% Lagrangian 
%%%function 
L = K_total - P_total; 

% dL/d(q_dot)
dL_dxdot = diff(L, x_dot);
dL_dtheta1dot = diff(L, theta_1_dot);
dL_dtheta2dot = diff(L, theta_2_dot);

% dL/dq
dL_dx = diff(L,x);
dL_dtheta_1 = diff(L, theta_1);
dL_dtheta_2 = diff(L, theta_2);

% dV/dq_dot
dV_dxdot = diff(V_total,x_dot);
dV_dtheta1dot = diff(V_total, theta_1_dot);
dV_dtheta2dot = diff(V_total, theta_2_dot);

% Lagrange's equations of motion
difeq_1 = diff(dL_dxdot, t) - dL_dx + dV_dxdot;
difeq_2 = diff(dL_dtheta1dot, t) - dL_dtheta_1 + dV_dtheta1dot;
difeq_3 = diff(dL_dtheta2dot, t) - dL_dtheta_2 + dV_dtheta2dot;

%% Solve the system of differential equations 

deqs = struct('x', difeq_1, 'theta_1', difeq_2, 'theta_2', difeq_3); 
var = [x;theta_1;theta_2];

leqs = [ deqs.x ==0; deqs.theta_1 == 0; deqs.theta_2 == 0];
%reduce order of differential equation
[eqs,vars] = reduceDifferentialOrder(leqs,var);

[MassMatrix,ForceMatrix] = massMatrixForm(eqs,vars);
% % MassMatrix = simplify(MassMatrix);
% % ForceMatrix = simplify(ForceMatrix);

MM = odeFunction(MassMatrix, vars);
FF = odeFunction(ForceMatrix, vars);

time = linspace(0, 60, 2000);
% initial conditions [x, theta1, theta2, xdot, theta1dot, theta2dot]
x_0 = [10 -pi/6+pi/3 pi/6+pi/3  0 0 0]; 
opts=odeset('Mass', MM, 'Stats','on');
[~, q] = ode45(FF, time, x_0, opts); % q: generalized coordinates

% Obtain the positions in Cartesian coordinates from the generalized
% coordinates
x = q(:, 1);
theta_1 = q(:, 2);
theta_2 = q(:, 3);

%%% m1     
x1 = x - l_1*sin(theta_1);
y1 = -l_1*cos(theta_1);

%%% m2 
x2 = x + l_1*sin(theta_2);
y2 = -l_2*cos(theta_2); 

%%% m3
d1 = l_1*cos((theta_2 - theta_1)/2);
A = (sin((theta_2 - theta_1)/2));
d2 = sqrt(l_2^2 - l_1^2*A.^2);
L = d1 + d2;
theta3 = (theta_1 + theta_2)/2;
S =sin(theta3);
C = cos(theta3);
x3 = x + L.*S;
y3 = -L.*C;
as = (d1 + d2/2);
ab = sqrt( (l_2/2)^2 - d2.^2*(1/4) + as.^2);
%% Simulate results to validate model
set(gcf, 'color', 'w')
set(gcf, 'position', [10, 100, 750, 750])

h = plot([]);
hold on
box on 
axis equal

for i = 1 : numel(time)
    if ~ishghandle(h)
        break
    end
    cla
    plot([0, x(i)], [0, 0], 'k', 'Linewidth', 2);
    plot(x(i), 0, '-', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
    plot([-x1(1), x1(i)], [0, y1(i)], 'k', 'Linewidth', 2);
    plot(x1(i), y1(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', m );
    plot([x2(1), x2(i)], [0, y2(i)], 'k', 'Linewidth', 2);
    plot(x2(i), y2(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', m );
    plot([0, x3(i)], [-ab(1), y3(i)], 'k', 'Linewidth', 2);
    plot(x3(i), y3(i), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', m );
    axis([-12, 12, -10, 0]);
    h = draw_spring_2D([x1(i); y1(i)], [x2(i); y2(i)], 10, 0.5);
    drawnow
end

function h = draw_spring_2D(A, B, number_of_coils, y_amplitude)    
    persistent t
    
    normalvector_AB = (B - A) / norm(B - A);
    offset_A = A + 1.25 * normalvector_AB;
    offset_B = B - 1.25 * normalvector_AB;
    distance_between_offsets = norm(offset_B - offset_A);
    
    t = linspace(-pi, number_of_coils * 2 * pi, 500);
    x_coordinate_between_offsets = distance_between_offsets * linspace(0, 1, numel(t));
    
    % ratio between x amplitude and y
    ratio_X_div_Y = 0.5;
    
    x = x_coordinate_between_offsets + ratio_X_div_Y * y_amplitude * cos(t);
    y = y_amplitude * sin(t);
    
    coil_positions = [x; y];
    rotation_matrix = [normalvector_AB, null(normalvector_AB')];
    rotated_coil_positions = rotation_matrix * coil_positions;
    h = plot([A(1), offset_A(1) + rotated_coil_positions(1,:), B(1)], ...
         [A(2), offset_A(2) + rotated_coil_positions(2,:), B(2)], 'k');
end

%% INPUTS
function [v, l_1, l_2, Mc, m, M1, M2, J1, J2, g, ks, b1, Fbr, bj ] = inputs
%INPUTS This function gives back the user defined parameters of the system 
%   The user-defined parameters are the following :
%   - DC motor voltage
%   - Rods' Lengths
%   - Masses
%   - Rotational Inertia
%   - Gravity acceleration
%   - Spring's stiffness
%   - Mechanical losses
%     * Ground friction
%     * Hydraulic brake force
%     * Joint friction
%   - Electrical losses
%     * At the DC motor

%% System parameters 

%%% Voltage
v = 15;%V

%%% Rods' Lengths l_1 = l_2, l_3 = l_4
l_1 = .3;%m %length of links 1,2
l_2 = .2;%m %length of links 3,4

%%%Masses
% Cart Mass Mc
Mc = 2;%kg

% Point masses m_1 = m_2 = m_3 = m
m = .5;%kg

%Rods' masses M and rotational inertia J
M1 = 0.12;%kg
M2 = 0.08;%kg
J1 = (1/3)*M1*(l_1/2)^2; %rotational inertia of links 1,2
J2 = (1/3)*M2*(l_2/2)^2; %rotational inertia of links 3,4

%%% Gravity acceleration 
g = 9.81; % m/s^2

%%% Spring's stiffness
ks = 50;

%%% Mechanical losses
% Ground friction
b1 = .2; % b1: friction coefficient
% Hydraulic brake force
Fbr = 200; %N
% At the joints
bj = 0.001;%N


%%% Electrical losses
% at the motor???
end