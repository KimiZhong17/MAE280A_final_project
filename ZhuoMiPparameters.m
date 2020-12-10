% best guess of Zhu Zhuo's calculation of MiP state-space model parameters
% with thanks to Matthew Pearson
% November 26, 2019

% physical constants
r = 0.034; % wheel radius m
l = 0.036; % center of mass to axle m
m_w = 0.027; % wheel mass kg
m_b = 0.263; % body mass kg
G_r = 35.57; % gearbox ratio
I_m = 3.6e-8; % motor armature moment of inertia kgm^2
I_b = 4e-4; % body moment of inertia
V_max = 7.4; % max battery voltage V
s = 0.003; % motor stall torque Nm
omega_f = 1760; % motor free run speed rad/sec
g = 9.8; % acceleration due to gravity m/s^2

% derived quantities
I_w = m_w*r^2/2+G_r^2*I_m; % momnet f interia of singel wheel plus gearbox
k = s/omega_f; % motor constant

% intermediate cosntants
a = 2*I_w+(m_b+2*m_w)*r^2;
b = m_b*r*l;
c = I_b+m_b*l^2;
d = m_b*g*l;
e = 2*G_r*s/V_max;
j = 2*G_r^2*k;
XX = 1/(a*c-b*b);

% Linearized matrices
A = [-(a+b)*j*XX (a+b)*j*XX a*d*XX;(b+c)*j*XX -(b+c)*j*XX -b*d*XX;1 0 0]
B = [-(a+b)*e*XX; (b+c)*e*XX; 0]
C = [1 0 0;0 1 0]
D = [0; 0]

% Zhu Zhuo's second-order discrete controller
ALCBKDmIp = [0.9129 0.0378;-0.0655 0.9933]
LLDmIp = [0.002835 -0.0005967;0.001463 -0.001286]
KLDmIp = [-380.821 322.4448]
DLDmIp = [-1.24321 0]