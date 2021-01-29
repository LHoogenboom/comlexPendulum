%Constants 
g = 9.81; %Gravity
rho = 1000; % Density

% Geometry 

L1 = 0.3; %Top Bar Length
L2 = 0.2; %Bottom Bar Length
r = 0.05; % Radius of the wheel 

% Masses

m0 = 2; %Cart
m1 = 0.5; %Mass 1
m2 = 0.5; %Mass 2
m3 = 0.5; %Mass 3
mL1 = L1*0.02*0.02*rho; % Top Bar Masses
mL2 = L2*0.02*0.02*rho; % Bottom Bar Masses
m_w = 0.5; % Mass of the wheel
J_w = 1/2*m_w*r^2; % Inertia of the wheel


% Spring 

d0 = 0.3; % rest length spring
k = 50; % spring constant

% Damping coefficients

b1 = 1e-3; %Damping of the links
b2 = 10; %Damping of the cart
b_w = 1e-3; %Friction of the wheels

% Motor 

L = 2.7e-3; %Indcutance
Kr = 0.7; %30.3e-3; %Motor Constant
R = 0.54; % Resistance

% Brake 
mus = 0.6; %Static friction coefficient
mudyn = 0.45; %Dynamic Friction coefficient
A1 = 2.5e-3; % Surface of the brake bad side
A2 = pi*(0.03/2)^2; %Surface of the disc brake side
F_pedal = 200; %Pedal Force
