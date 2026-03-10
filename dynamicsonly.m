%% dynamics_only_yokota_hcw.m
% Dynamics-only model based on:
% Yokota et al., "Differential Drag Estimation and Autonomous Control
% for Non-Cooperative Rendezvous"
%
% Paper-consistent circular-orbit relative motion:
%   xdd =  2*w*zd
%   ydd = -w^2*y
%   zdd = -2*w*xd + 3*w^2*z
%
% State vector:
%   X = [x; y; z; xd; yd; zd]

clear; clc; close all;

%% 1) Orbit parameters
mu = 3.986004418e14;      % [m^3/s^2] Earth gravitational parameter
Re = 6378.137e3;          % [m] Earth radius
alt = 400e3;              % [m] reference orbit altitude
r0 = Re + alt;            % [m] orbital radius

w = sqrt(mu / r0^3);      % [rad/s] mean motion
Torb = 2*pi / w;          % [s] orbital period

fprintf('Mean motion w  = %.6e rad/s\n', w);
fprintf('Orbital period = %.2f min\n', Torb/60);

%% 2) Paper-consistent continuous-time dynamics
% Xdot = A*X + B*u
% For dynamics-only, u = 0

A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 0 0 2*w;
     0 -w^2 0 0 0 0;
     0 0 3*w^2 -2*w 0 0];

B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];

%% 3) Simulation settings
nOrbits = 4;
dt = 0.5;                         % [s]
t = (0:dt:nOrbits*Torb).';
N = numel(t);

% Exact discrete-time propagation for zero-order hold
M = expm([A B; zeros(3,9)] * dt);
Ad = M(1:6,1:6);

%% 4) Initial condition
% Keep this small and simple for baseline validation
% X0 = [x0; y0; z0; xd0; yd0; zd0]

X = zeros(6, N);
X(:,1) = [0; 0; -30; 0; 0; 0];

%% 5) Propagate with zero input
for k = 1:N-1
    X(:,k+1) = Ad * X(:,k);
end

%% 6) Extract states
x  = X(1,:);
y  = X(2,:);
z  = X(3,:);
xd = X(4,:);
yd = X(5,:);
zd = X(6,:);

%% 7) Plot time histories
figure();
plot(t/Torb, x, 'LineWidth', 1.4); hold on;
plot(t/Torb, y, 'LineWidth', 1.4);
plot(t/Torb, z, 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('Position [m]');
title('Dynamics-only relative motion: positions');
legend('x','y','z','Location','best');

figure();
plot(t/Torb, xd, 'LineWidth', 1.4); hold on;
plot(t/Torb, yd, 'LineWidth', 1.4);
plot(t/Torb, zd, 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('Velocity [m/s]');
title('Dynamics-only relative motion: velocities');
legend('xd','yd','zd','Location','best');

%% 8) Plot geometry
figure();
plot3(y, z, x, 'LineWidth', 1.5);
grid on;
xlabel('y [m]');
ylabel('z [m]');
zlabel('x [m]');
title('Relative trajectory in Hill coordinates');
view(35,20);

figure();
plot(z, x, 'LineWidth', 1.5);
grid on;
xlabel('z [m]');
ylabel('x [m]');
title('x-z plane trajectory');

figure();
plot(y, x, 'LineWidth', 1.5);
grid on;
xlabel('y [m]');
ylabel('x [m]');
title('x-y plane trajectory');

%% 9) Optional analytical consistency check
% For dynamics-only with zero input, X(t) = Phi(t)*X0
% Compare numerical propagation vs analytic Phi at final time.

X0 = X(:,1);
Phi_final = Phi_paper(t(end), w);
X_final_analytic = Phi_final * X0;
X_final_numeric  = X(:,end);

disp('Final state error (numeric - analytic):');
disp(X_final_numeric - X_final_analytic);

%% Local function: Phi(t) from paper Eq. (16)
function Phi = Phi_paper(tau, w)
    c = cos(w*tau);
    s = sin(w*tau);

    Phi = [ ...
        1, 0, 6*(w*tau - s), (4*s - 3*w*tau)/w, 0, 2*(1-c)/w;
        0, c, 0,            0,                 s/w, 0;
        0, 0, 4 - 3*c,      0,                -2*(1-c)/w, s/w;
        0, 0, 6*w*(1-c),    4*c - 3,           0, 2*s;
        0, -w*s, 0,         0,                 c, 0;
        0, 0, 3*w*s,       -2*s,               0, c];
end