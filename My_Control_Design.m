%% MAE280A Final Project
% Written by Guoren Zhong
% gzhong@eng.ucsd.edu
% Dec.10, 2020
clc;
clear all;

%% Parameters
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
A = [-(a+b)*j*XX (a+b)*j*XX a*d*XX;(b+c)*j*XX -(b+c)*j*XX -b*d*XX;1 0 0];
B = [-(a+b)*e*XX; (b+c)*e*XX; 0];
C = [1 0 0;0 1 0];
D = [0; 0];

% Zhu Zhuo's second-order discrete controller
ALCBKDmIp = [0.9129 0.0378;-0.0655 0.9933];
LLDmIp = [0.002835 -0.0005967;0.001463 -0.001286];
KLDmIp = [-380.821 322.4448];
DLDmIp = [-1.24321 0];

i = sqrt(-1);

%% Basic design 
% capability: thetaic_max = 0.36
%             (single) thetadot_noise_max = 0.0094 (in 30 sec)
%             (single) phidot_noise_max = 0.024 (in 30 sec)
%             thetaic_max = 0.67 (with the initial condition for controller
%                                 being [5.3 0 0])
%             thetaic_max = 0.65 (with the initial condition for thetadot
%                                 being -5)
mipC = ss(A,B,C,D);                 % continuous-time system
mipD = c2d(mipC,0.01);              % convert to discrete-time system with 100Hz
[AD, BD, CD, DD] = ssdata(mipD);
Kmine = place(AD,BD,[0.95;0.94;0.93]);
Lmine = place(AD',CD',[0.95;0.94;0.93])';
Amine = AD - BD * Kmine - Lmine * CD;
Dmine = [0 0];
thetaic = 0.36;
%% Get data for controller 1 and plot
t1 = out.theta.Time;
thetadot1 = out.ScopeData{1}.Values.Data;
phidot1 = out.ScopeData{2}.Values.Data;
theta1 = out.ScopeData{3}.Values.Data;
thetadot1_est = out.ScopeData{4}.Values.Data;
phidot1_est = out.ScopeData{5}.Values.Data;
theta1_est = out.ScopeData{6}.Values.Data;
u1_computed = out.u_computed.Data;
u1_applied = out.u_applied.Data;

figure(1)
plot(t1,thetadot1,':',t1,phidot1,':',t1,theta1,':','LineWidth',1.5)
grid on
grid minor
xlabel('Time');
ylabel('rad & rad/s')
set(gcf,'unit','centimeters','position',[6,10,15,6])
set(gca,'position',[0.09,0.18,0.89,0.79])
hold on
stairs(t1,[thetadot1_est,phidot1_est,theta1_est],'LineWidth',1.5)
leg1 = legend('$\dot{\theta}$ (rad/s)','$\dot{\phi}$ (rad/s)','$\theta$ (rad)',...
    '$\dot{\theta}_{est}$ (rad/s)','$\dot{\phi}_{est}$ (rad/s)','$\theta_{est}$ (rad)');
set(leg1, 'Interpreter', 'latex');
hold off

figure(2)
plot(t1(1:1001),u1_computed,'-',t1(1:1001),u1_applied,':','LineWidth',1.5)
grid on
grid minor
xlabel('Time');
ylabel('Voltage')
set(gcf,'unit','centimeters','position',[22,10,15,6])
set(gca,'position',[0.09,0.18,0.89,0.79])
legend('V_{computed}', 'V_{applied}')
%% Choose faster poles 
% capability: thetaic_max = 0.57
%             (single) thetadot_noise_max = 0.035 (in 30 sec)
%             (single) phidot_noise_max = 0.065 (in 30 sec)
%             (combine) noise_max = 0.092
%             thetaic_max = 1.09 (with the initial estimation for thetadot
%                                 being 9.7)
%             thetaic_max = 0.76 (with the initial condition for thetadot
%                                 being -8.5)
%             thetaic_max = 1.03 (with the initial estimation for theta
%                                 being 0.55)
mipC = ss(A,B,C,D);                 % continuous-time system
mipD = c2d(mipC,0.01);              % convert to discrete-time system with 100Hz
[AD, BD, CD, DD] = ssdata(mipD);
Kmine = place(AD,BD,[0.9;0.88;0.86]);       % state feedback
Lmine = place(AD',CD',[0.9;0.88;0.86])';    % controller
Amine = AD - BD * Kmine - Lmine * CD;
Dmine = [0 0];
thetaic = 1.03;
%% Get data for controller 2 and plot
t2 = out.theta.Time;
thetadot2 = out.ScopeData{1}.Values.Data;
phidot2 = out.ScopeData{2}.Values.Data;
theta2 = out.ScopeData{3}.Values.Data;
thetadot2_est = out.ScopeData{4}.Values.Data;
phidot2_est = out.ScopeData{5}.Values.Data;
theta2_est = out.ScopeData{6}.Values.Data;
u2_computed = out.u_computed.Data;
u2_applied = out.u_applied.Data;

figure(1)
plot(t2,thetadot2,':',t2,phidot2,':',t2,theta2,':','LineWidth',1.5)
grid on
grid minor
xlabel('Time');
ylabel('rad & rad/s')
%ylim([-2,4.5])
set(gcf,'unit','centimeters','position',[6,10,15,6])
set(gca,'position',[0.09,0.18,0.89,0.79])
hold on
stairs(t2,[thetadot2_est,phidot2_est,theta2_est],'LineWidth',1.5)
leg2 = legend('$\dot{\theta}$ (rad/s)','$\dot{\phi}$ (rad/s)','$\theta$ (rad)',...
    '$\dot{\theta}_{est}$ (rad/s)','$\dot{\phi}_{est}$ (rad/s)','$\theta_{est}$ (rad)');
set(leg2, 'Interpreter', 'latex');
hold off

figure(2)
plot(t2(1:201),u2_computed,'-',t2(1:201),u2_applied,':','LineWidth',1.5)
grid on
grid minor
xlabel('Time');
ylabel('Voltage')
%ylim([4,7])
%xlim([0.05,0.2])
set(gcf,'unit','centimeters','position',[22,10,15,6])
set(gca,'position',[0.09,0.18,0.89,0.79])
legend('V_{computed}', 'V_{applied}')

%% Compute K by solving the Discrete Algebraic Riccati Equation
% capability: thetaic_max = 0.65
%             noise_power_max = 0.031
%             (single) thetadot_noise_max = 0.04 (in 30 sec)
%             (single) phidot_noise_max = 0.04 (in 30 sec)
% improve by changing initial conditions:
%             thetaic_max = 0.95 (with the initial estimation for thetadot
%                                 being 3)
%             thetaic_max = 0.80 (with the initial condition for thetadot
%                                 being -6.5)
%             thetaic_max = 0.978 (with the initial estimation for theta
%                                 being 0.26)
mipC = ss(A,B,C,D);                 % continuous-time system
mipD = c2d(mipC,0.01);              % convert to discrete-time system with 100Hz
[AD, BD, CD, DD] = ssdata(mipD);
[~,~,Kmine] = dare(AD,BD,diag([0 1 0]),0.01);   % state feedback
Lmine = place(AD',CD',[0.9;0.88;0.86])';        % controller
Amine = AD - BD * Kmine - Lmine * CD;
Dmine = [0 0];
thetaic = 0.978;
%% Get data for controller 3 and plot
t3 = out.theta.Time;
thetadot3 = out.ScopeData{1}.Values.Data;
phidot3 = out.ScopeData{2}.Values.Data;
theta3 = out.ScopeData{3}.Values.Data;
thetadot3_est = out.ScopeData{4}.Values.Data;
phidot3_est = out.ScopeData{5}.Values.Data;
theta3_est = out.ScopeData{6}.Values.Data;
u3_computed = out.u_computed.Data;
u3_applied = out.u_applied.Data;

figure(1)
plot(t3,thetadot3,':',t3,phidot3,':',t3,theta3,':','LineWidth',1.5)
grid on
grid minor
xlabel('Time');
ylabel('rad & rad/s')
%ylim([-2,4.5])
set(gcf,'unit','centimeters','position',[6,10,15,6])
set(gca,'position',[0.09,0.18,0.89,0.79])
hold on
stairs(t3,[thetadot3_est,phidot3_est,theta3_est],'LineWidth',1.5)
leg3 = legend('$\dot{\theta}$ (rad/s)','$\dot{\phi}$ (rad/s)','$\theta$ (rad)',...
    '$\dot{\theta}_{est}$ (rad/s)','$\dot{\phi}_{est}$ (rad/s)','$\theta_{est}$ (rad)');
set(leg3, 'Interpreter', 'latex');
hold off

figure(2)
plot(t3(1:201),u3_computed,'-',t3(1:201),u3_applied,':','LineWidth',1.5)
grid on
grid minor
xlabel('Time');
ylabel('Voltage')
%xlim([0,0.15])
%ylim([5,7])
set(gcf,'unit','centimeters','position',[22,10,15,6])
set(gca,'position',[0.09,0.18,0.89,0.79])
legend('V_{computed}', 'V_{applied}')