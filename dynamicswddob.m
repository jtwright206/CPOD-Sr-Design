%% ddob_rls_yokota_hcw.m
% Dynamics + isolated DDOB estimator block + RLS atmospheric density estimator
% Based on:
% Yokota et al., "Differential Drag Estimation and Autonomous Control
% for Non-Cooperative Rendezvous"
%
% This script includes:
%   1) Paper-consistent HCW/Hill relative motion
%   2) A truth differential-drag disturbance input
%   3) An isolated DDOB estimator block
%   4) An RLS atmospheric-density estimator
%
% It does NOT yet include:
%   - guidance
%   - control
%   - feedback trigger

clear; clc; close all;

%% 1) Orbit parameters
mu = 3.986004418e14;      % [m^3/s^2]
Re = 6378.137e3;          % [m]
alt = 400e3;              % [m]
r0 = Re + alt;            % [m]

w = sqrt(mu / r0^3);      % [rad/s]
Torb = 2*pi / w;          % [s]
Vt = sqrt(mu / r0);       % [m/s] target orbital speed magnitude

fprintf('Mean motion w  = %.6e rad/s\n', w);
fprintf('Orbital period = %.2f min\n', Torb/60);
fprintf('Orbital speed  = %.3f km/s\n', Vt/1000);

%% 2) Spacecraft / drag parameters (for paper-style drag model)
Cd = 2.5;
mc = 500;     % [kg] chaser mass
Ac = 1;       % [m^2] chaser area
mt = 3000;    % [kg] target mass
At = 40;      % [m^2] target area

% Useful scalar coefficient in the paper's drag model
% d = rho * phi
% where phi depends on relative speed and spacecraft parameters
% (Eqs. 20-26 concept)
% We'll compute phi(k) each step.

%% 3) Paper-consistent continuous-time dynamics
% State X = [x; y; z; xd; yd; zd]
%
% xdd =  2*w*zd + ax
% ydd = -w^2*y  + ay
% zdd = -2*w*xd + 3*w^2*z + az

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

%% 4) Simulation settings
nOrbits = 4;
dt = 0.1;                         % [s]
t = (0:dt:nOrbits*Torb).';
N = numel(t);

% Exact discrete-time propagation
M = expm([A B; zeros(3,9)] * dt);
Ad = M(1:6,1:6);
Bd = M(1:6,7:9);

%% 5) Initial condition
X = zeros(6, N);
X(:,1) = [0; 0; -30; 0; 0; 0];

%% 6) Truth density model
% Instead of prescribing drag magnitude directly, prescribe density and
% generate drag through the paper-style drag model.
rhoTrue = zeros(1, N);

for k = 1:N
    tk = t(k);

    % Smooth varying density truth
    rhoTrue(k) = 4e-12 * (1 + 0.3*sin(2*pi*tk/Torb));   % [kg/m^3]
end

%% 7) Truth disturbance model from density
% dTruth(:,k) = rhoTrue(k) * phiTruth(k) * drag_direction
%
% Keep drag along x only for now, as in your current simplified setup.

dTruth   = zeros(3, N);
phiTruth = zeros(1, N);

for k = 1:N
    vrel_mag = norm(X(4:6, max(k,1)));  % initially zero, later overwritten more meaningfully
    phiTruth(k) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag)^2 ...
                - 0.5*Cd*(At/mt)*(Vt)^2;

    dTruth(:,k) = [rhoTrue(k) * phiTruth(k);
                   0;
                   0];
end

%% 8) Simulate truth dynamics with zero control
u = zeros(3,1);

for k = 1:N-1
    % recompute phi from current truth state velocity
    vrel_mag = norm(X(4:6,k));
    phiTruth(k) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag)^2 ...
                - 0.5*Cd*(At/mt)*(Vt)^2;

    dTruth(:,k) = [rhoTrue(k) * phiTruth(k);
                   0;
                   0];

    X(:,k+1) = Ad*X(:,k) + Bd*(u + dTruth(:,k));
end

% fill final sample consistently
vrel_mag = norm(X(4:6,N));
phiTruth(N) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag)^2 ...
            - 0.5*Cd*(At/mt)*(Vt)^2;
dTruth(:,N) = [rhoTrue(N) * phiTruth(N); 0; 0];

%% 9) Measurement model
% Start with noise-free state for DDOB + RLS validation.

Xmeas = X;

%% 10) Build measured acceleration from measured velocity
accMeas = zeros(3, N);

for k = 2:N
    accMeas(:,k) = (Xmeas(4:6,k) - Xmeas(4:6,k-1)) / dt;
end

%% 11) Isolated DDOB block (paper Eq. 18 concept)
% DDOB raw estimate:
%   d_hat_raw = xdd_meas - model_term - u - a_d
%
% Here:
%   u = 0
%   a_d = 0
%
% Then apply low-pass filter Q.

fc_Q = 10;                          % [Hz], paper Table 1 value
alphaQ = exp(-2*pi*fc_Q*dt);

dHatRaw  = zeros(3, N);
dHatFilt = zeros(3, N);

for k = 1:N
    xm = Xmeas(:,k);

    % model term from paper's circular HCW dynamics
    modelTerm = [ 2*w*xm(6);
                 -w^2*xm(2);
                 -2*w*xm(4) + 3*w^2*xm(3) ];

    dHatRaw(:,k) = accMeas(:,k) - modelTerm - u;

    if k == 1
        dHatFilt(:,k) = dHatRaw(:,k);
    else
        dHatFilt(:,k) = alphaQ*dHatFilt(:,k-1) + (1-alphaQ)*dHatRaw(:,k);
    end
end

%% 12) RLS atmospheric-density estimator
% Paper structure:
%   y = |DDOB|
%   y = phi * rho
% Estimate rho with scalar RLS, then reconstruct d_hat = rho_hat * phi.
%
% Since our simplified truth drag is along x only, use x-component here.

lambdaRLS = 0.999;
P = 1e20;                 % large initial covariance
rhoHat = zeros(1, N);
rhoHat(1) = 3e-12;        % initial guess

phiMeas = zeros(1, N);
phiFilt = zeros(1, N);
yRLS    = zeros(1, N);
dHatRLS = zeros(3, N);

for k = 1:N
    % Measured regressor phi from measured relative speed
    vrel_mag_meas = norm(Xmeas(4:6,k));
    phiMeas(k) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag_meas)^2 ...
               - 0.5*Cd*(At/mt)*(Vt)^2;

    % Paper uses Q in the regression path for phase alignment
    if k == 1
        phiFilt(k) = phiMeas(k);
    else
        phiFilt(k) = alphaQ*phiFilt(k-1) + (1-alphaQ)*phiMeas(k);
    end

    % Measurement y = |filtered DDOB|
    % For this simplified model, use x-component magnitude
    yRLS(k) = abs(dHatFilt(1,k));

    if k > 1
        phi_k = phiFilt(k);

        if abs(phi_k) > 1e-20
            K = (P * phi_k) / (lambdaRLS + phi_k^2 * P);
            rhoHat(k) = rhoHat(k-1) + K * (yRLS(k) - phi_k * rhoHat(k-1));
            P = (1/lambdaRLS) * (P - K * phi_k * P);
        else
            rhoHat(k) = rhoHat(k-1);
        end
    end

    % guard against nonphysical negative density
    rhoHat(k) = max(rhoHat(k), 0);

    % Reconstructed disturbance estimate from rho_hat
    dHatRLS(:,k) = [rhoHat(k) * phiMeas(k);
                    0;
                    0];
end

%% 13) Plot relative motion
figure();
plot3(X(2,:), X(3,:), X(1,:), 'LineWidth', 1.4);
grid on;
xlabel('y [m]');
ylabel('z [m]');
zlabel('x [m]');
title('Relative trajectory with truth disturbance');
view(35,20);

figure();
plot(X(3,:), X(1,:), 'LineWidth', 1.5);
grid on;
xlabel('z [m]');
ylabel('x [m]');
title('x-z trajectory with truth disturbance');

%% 14) Plot drag estimation comparison
figure();
plot(t/Torb, dTruth(1,:), 'k', 'LineWidth', 1.6); hold on;
plot(t/Torb, dHatFilt(1,:), 'r', 'LineWidth', 1.0);
plot(t/Torb, dHatRLS(1,:), 'b', 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('d_x [m/s^2]');
title('Differential drag estimation');
legend('Truth', 'DDOB filtered', 'DDOB + RLS', 'Location', 'best');

%% 15) Plot atmospheric density estimate
figure();
plot(t/Torb, rhoTrue, 'k', 'LineWidth', 1.6); hold on;
plot(t/Torb, rhoHat, 'b', 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('\rho [kg/m^3]');
title('Atmospheric density estimation (RLS)');
legend('Truth', '\rho-hat', 'Location', 'best');

%% 16) Simple estimation error summary
errRaw  = dHatRaw - dTruth;
errFilt = dHatFilt - dTruth;
errRLS  = dHatRLS - dTruth;

rmseRaw  = sqrt(mean(vecnorm(errRaw,2,1).^2));
rmseFilt = sqrt(mean(vecnorm(errFilt,2,1).^2));
rmseRLS  = sqrt(mean(vecnorm(errRLS,2,1).^2));

rhoRMSE = sqrt(mean((rhoHat - rhoTrue).^2));

fprintf('\nDDOB raw RMSE      = %.8e m/s^2\n', rmseRaw);
fprintf('DDOB filt RMSE     = %.8e m/s^2\n', rmseFilt);
fprintf('DDOB + RLS RMSE    = %.8e m/s^2\n', rmseRLS);
fprintf('Density RMSE       = %.8e kg/m^3\n', rhoRMSE);