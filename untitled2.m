clear; clc;
%% Inputs
[v, l_1, l_2, Mc, m, M1, M2, J1, J2, g, ks, b1, Fbr, bj] = inputs;

car = struct('V', v, 'm', Mc, 'x', 0, 'xd', 0, 'mu', b1);
spr = struct('k',ks);
p1 = struct('m', m, 'mu', bj);
p2 = struct('m', m, 'mu', bj);
p3 = struct('m', m, 'mu', bj);
l1 = struct('l',l_1,'m',M1,'J',J1);
l2 = l1;
l3 = struct('l',l_2,'m',M2,'J',J2);
l4 = l3;

%%% Generalized coordinates
% - x: position of the cart 
% - theta_1: rotational angle of rod 1 with respect to the vertical
% - theta_2: rotational angle of rod 2, with respect ot the vertical
syms x(t) theta_1(t)  theta_2(t)
car.x = str2sym('x(t)');
th1 = str2sym('theta_1(t)');
th2 = str2sym('theta_2(t)');
th3 = (th1+th2)/2;
th4 = (th2-th1)/2;

% kite height (L)
d1 = l1.l*cos(th4); % top iscoleces height
d2 = sqrt(l2.l^2 - l1.l^2*sin(th4)^2); % bot iscoleces height
L = d1 + d2; % Kite height (L)

clear l_1 l_2 v Mc m M1 M2 J1 J2 ks
%% Point Mass Coordinates
% Top of kite has no mass and always at (0,0)

%Point masses (p)
p1.x = car.x + l1.l*sin(th1);
p1.y = -l1.l*cos(th1);

p2.x = car.x + l1.l*sin(th2);
p2.y = -l2.l*cos(th2);

p3.x = car.x + L*sin(th3); 
p3.y = -L*cos(th3);

%% Rod Coordinates

% top left link
l1.x = car.x + (l1.l/2)*sin(th1);
l1.y = -(l1.l/2)*cos(th1);

% top right link
l2.x = car.x + (l2.l/2)*sin(th2);
l2.y = -(l2.l/2)*cos(th2);

% bot left link
BD = l3.l/2; % half link length
DN = d2/2; % follows from similar triangles (Angle-Angle theorem)
AN = L-DN;
BN = sqrt(BD^2-DN^2); % pyhtagoras
AB = sqrt(BN^2+AN^2);

% these lines are problematic
num = l1.l^2 + (l3.l/2)^2 - (sqrt(l3.l^2 - l1.l^2*(sin((th2 - th1)/2))^2)/2)^2 + (l1.l*cos((th2 - th1)/2) + sqrt(l3.l^2 - l1.l^2*(sin((th2 - th1)/2))^2)/2)^2 - (l3.l/2)^2;
den = 2*l1.l * sqrt((l3.l/2)^2 - (sqrt(l3.l^2 - l1.l^2*(sin((th2 - th1)/2))^2)/2)^2 + (l1.l*cos((th2 - th1)/2) + sqrt(l3.l^2 - l1.l^2*(sin((th2 - th1)/2))^2)/2)^2);
beta2 = acos(num/den);

th5 = th1 + beta2;

l3.x = car.x + AB*sin(th5);
l3.y = -AB*cos(th5);

% bot right link
AC = AB;
th6 = th2 - beta2;
l4.x = car.x + AC*sin(th6); 
l4.y = -AC*cos(th6);

clear num den
%% Lineal velocities of points & Links
car.xd = diff(car.x,t);
p1.xd = diff(p1.x,t); 
p1.yd = diff(p1.y,t);
p1.v = norm(p1.xd,p1.yd);

p2.xd = diff(p2.x,t); 
p2.yd = diff(p2.y,t);
p2.v = norm(p2.xd,p2.yd);

p3.xd = diff(p3.x,t); 
p3.yd = diff(p3.y,t);
p3.v = norm(p3.xd,p3.yd);

l1.xd = diff(l1.x,t);
l1.yd = diff(l1.y,t);
l1.v = norm(l1.xd,l1.y);

l2.xd = diff(l2.x,t);
l2.yd = diff(l2.y,t);
l2.v = norm(l2.xd,l2.y);

l3.xd = diff(l3.x,t);
l3.yd = diff(l3.y,t);
l3.v = norm(l3.xd,l3.y);

l4.xd = diff(l4.x,t);
l4.yd = diff(l4.y,t);
l4.v = norm(l4.xd,l4.y);

%% Angular velocities of Links

th1d = diff(th1,t);
th2d = diff(th2,t);

eta1 = pi - acos((l2.l^2 + l1.l^2 - (l1.l*cos((th2 - th1)/2) + sqrt(l2.l^2 - l1.l^2*(sin((th2 - th1)/2))^2))^2)/(2*l2.l*l1.l));
eta1d = diff(eta1,t); 

eta2 = eta1; 
eta2d = diff(eta2,t);

%% Kinetic energies

%car
car.K = (1/2)*car.m*car.xd^2;

% Point masses
p1.K = (1/2)*p1.m*p1.v^2;
p2.K = (1/2)*p2.m*p2.v^2;
p3.K = (1/2)*p3.m*p3.v^2;

% Links (Lineal + Angular)
l1.K = (1/2)*l1.m*l1.v^2 + (1/2)*l1.J*th1d^2;
l2.K = (1/2)*l2.m*l2.v^2 + (1/2)*l2.J*th2d^2;
l3.K = (1/2)*l3.m*l3.v^2 + (1/2)*l3.J*eta1d^2;
l4.K = (1/2)*l4.m*l4.v^2 + (1/2)*l4.J*eta2d^2;

Ktot = car.K + p1.K + p2.K + p3.K + l1.K + l2.K + l3.K + l4.K;

%% Potential Energies

% car
car.P = 0;

% point masses
p1.P = p1.m*g*p1.y;
p2.P = p2.m*g*p2.y;
p3.P = p3.m*g*p3.y;

% Links
l1.P = l1.m*g*l1.y;
l2.P = l2.m*g*l2.y;
l3.P = l3.m*g*l3.y;
l4.P = l4.m*g*l4.y;

% spring
spr.l = abs(norm(p2.x-p1.x , p2.y-p1.y)); % abs length (l0 = 0)
spr.P = (1/2)*spr.k*spr.l^2;

Ptot = car.P + p1.P + p2.P + p3.P + l1.P + l2.P + l3.P + l4.P + spr.P; 

%% Dissipatives energy
% Assuming the system ais always moving, a constant dynamic friction
% coefficient can be used.
car.D = car.m*car.mu*g;
p1.D = p1.m*p1.mu*g;
p2.D = p2.m*p2.mu*g;
p3.D = p3.m*p3.mu*g;

Dtot  = car.D + p1.D + p2.D + p3.D;

%% Lagrangian 
%%%function 
L = Ktot - Ptot; 

% dL/d(q_dot)
dL_dxdot = diff(L, car.xd);
dL_dtheta1dot = diff(L, th1d);
dL_dtheta2dot = diff(L, th2d);

% dL/dq
dL_dx = diff(L,x);
dL_dtheta_1 = diff(L, th1);
dL_dtheta_2 = diff(L, th2);

% dV/dq_dot
dV_dxdot = diff(Dtot,car.xd);
dV_dtheta1dot = diff(Dtot, th1d);
dV_dtheta2dot = diff(Dtot, th2d);

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
x1 = car.x - l1.l*sin(theta_1);
y1 = -l1.l*cos(theta_1);

%%% m2 
x2 = car.x + l1.l*sin(theta_2);
y2 = -l2.l*cos(theta_2); 

%%% m3
d1 = l1.l*cos((theta_2 - theta_1)/2);
A = (sin((theta_2 - theta_1)/2));
d2 = sqrt(l2.l^2 - l1.l^2*A.^2);
L = d1 + d2;
theta3 = (theta_1 + theta_2)/2;
S =sin(theta3);
C = cos(theta3);
x3 = car.x + L.*S;
y3 = -L.*C;
as = (d1 + d2/2);
ab = sqrt( (l2.l/2)^2 - d2.^2*(1/4) + as.^2);
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


%% Functions

function [v, l1, l2, Mc, m, M1, M2, J1, J2, g, ks, b1, Fbr, bj ] = inputs
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
l1 = 2;%m %length of links 1,2
l2 = 1.2;%m %length of links 3,4

%%%Masses
% Cart Mass Mc
Mc = 100;%kg

% Point masses m_1 = m_2 = m_3 = m
m = 10;%kg

%Rods' masses M and rotational inertia J
M1 = 1;%kg
M2 =0.8;%kg
J1 = (1/3)*M1*l1^2; %rotational inertia of links 1,2
J2 = (1/3)*M2*l2^2; %rotational inertia of links 3,4

%%% Gravity acceleration 
g = 10; % m/s^2

%%% Spring's stiffness
ks = 0.5;

%%% Mechanical losses
% Ground friction
b1 = 0.45; % b1: friction coefficient
% Hydraulic brake force
Fbr = 1000; %N
% At the joints
bj = 0.5;%N


%%% Electrical losses
% at the motor???
end