%% ddob_rls_guidance_control_yokota_hcw.m
% Dynamics + DDOB + RLS + Guidance + Control
% Based on:
% Yokota et al., "Differential Drag Estimation and Autonomous Control
% for Non-Cooperative Rendezvous"
%
% This script includes:
%   1) Paper-consistent HCW/Hill relative motion
%   2) Truth differential-drag disturbance input
%   3) DDOB estimator block
%   4) RLS atmospheric-density estimator
%   5) Guidance block based on Phi/Gamma and handover target mapping
%   6) PD tracking controller with saturation
%
% Still simplified relative to the full paper:
%   - no J2 modeled disturbance
%   - no optimization over ROE target bounds
%   - no FB ON/OFF trigger block yet
%   - drag kept along x only for clarity

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

%% 2) Spacecraft / drag parameters
Cd = 2.5;
mc = 500;     % [kg] chaser mass
Ac = 1;       % [m^2] chaser area
mt = 3000;    % [kg] target mass
At = 40;      % [m^2] target area

%% 3) Controller / estimator parameters
fc_Q = 10;                    % [Hz]
alphaQ = exp(-2*pi*fc_Q*dt_placeholder(0.1));  % temporary placeholder, overwritten later

lambdaRLS = 0.999;
Kp = 5e-4;                    % reduced from paper for stable first implementation
Kd = 5e-3;
Tmax = 8;                     % [N]
umax = Tmax / mc;             % [m/s^2]

%% 4) Paper-consistent continuous-time dynamics
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

%% 5) Simulation settings
nOrbits = 4;
dt = 0.1;                         % [s]
t = (0:dt:nOrbits*Torb).';
N = numel(t);

% update Q filter coefficient now that dt is known
alphaQ = exp(-2*pi*fc_Q*dt);

% Exact discrete-time propagation
M = expm([A B; zeros(3,9)] * dt);
Ad = M(1:6,1:6);
Bd = M(1:6,7:9);

%% 6) Guidance settings
% Paper Eq. (39) target vector:
% a * delta_alpha_target = [-10; -300; 0; -100; 0; -100]
%
% To keep motion nearly planar at first, you may later set the last term to 0.
a_delta_alpha_target = [-10; -300; 0; -100; 0; -100];

% Guidance segment length:
% shorter horizon prevents huge demanded corrections
tau_handover = 300;              % [s]

segment_active = false;
segment_t0 = 0;
segment_tau = tau_handover;
x_plan0 = zeros(6,1);
d_plan  = zeros(3,1);

%% 7) Initial condition
X = zeros(6, N);

% Start near the desired trajectory but with an offset
x_target_init = target_state_from_roe(a_delta_alpha_target, tau_handover, w);
X(:,1) = x_target_init + [-20; 0; -10; 0; 0; 0];

%% 8) Truth density model
rhoTrue = zeros(1, N);

for k = 1:N
    tk = t(k);
    rhoTrue(k) = 4e-12 * (1 + 0.3*sin(2*pi*tk/Torb));   % [kg/m^3]
end

%% 9) Allocate truth / estimator / controller arrays
dTruth   = zeros(3, N);
phiTruth = zeros(1, N);

Xmeas    = zeros(6, N);
accMeas  = zeros(3, N);

dHatRaw  = zeros(3, N);
dHatFilt = zeros(3, N);
dHatRLS  = zeros(3, N);

rhoHat   = zeros(1, N);
rhoHat(1) = 3e-12;
phiMeas  = zeros(1, N);
phiFilt  = zeros(1, N);
yRLS     = zeros(1, N);
P = 1e20;

xGuid    = zeros(6, N);
uCmd     = zeros(3, N);
trackErr = zeros(1, N);

u = zeros(3,1);

%% 10) Main simulation loop
for k = 1:N-1
    tk = t(k);

    % ------------------------------------------------------------
    % 10.1 Truth differential drag
    % ------------------------------------------------------------
    vrel_mag_truth = norm(X(4:6,k));
    phiTruth(k) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag_truth)^2 ...
                - 0.5*Cd*(At/mt)*(Vt)^2;

    dTruth(:,k) = [rhoTrue(k) * phiTruth(k);
                   0;
                   0];

    % ------------------------------------------------------------
    % 10.2 Measurement model
    % ------------------------------------------------------------
    % Keep noise-free for now
    Xmeas(:,k) = X(:,k);

    if k >= 2
        accMeas(:,k) = (Xmeas(4:6,k) - Xmeas(4:6,k-1)) / dt;
    end

    % ------------------------------------------------------------
    % 10.3 DDOB block
    % ------------------------------------------------------------
    xm = Xmeas(:,k);

    modelTerm = [ 2*w*xm(6);
                 -w^2*xm(2);
                 -2*w*xm(4) + 3*w^2*xm(3) ];

    dHatRaw(:,k) = accMeas(:,k) - modelTerm - u;

    if k == 1
        dHatFilt(:,k) = dHatRaw(:,k);
    else
        dHatFilt(:,k) = alphaQ*dHatFilt(:,k-1) + (1-alphaQ)*dHatRaw(:,k);
    end

    % ------------------------------------------------------------
    % 10.4 RLS atmospheric-density estimator
    % ------------------------------------------------------------
    vrel_mag_meas = norm(Xmeas(4:6,k));
    phiMeas(k) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag_meas)^2 ...
               - 0.5*Cd*(At/mt)*(Vt)^2;

    if k == 1
        phiFilt(k) = phiMeas(k);
    else
        phiFilt(k) = alphaQ*phiFilt(k-1) + (1-alphaQ)*phiMeas(k);
    end

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

    rhoHat(k) = max(rhoHat(k), 0);

    dHatRLS(:,k) = [rhoHat(k) * phiMeas(k);
                    0;
                    0];

    % ------------------------------------------------------------
    % 10.5 Guidance planning
    % ------------------------------------------------------------
    need_replan = false;

    if ~segment_active
        need_replan = true;
    elseif (tk - segment_t0) >= segment_tau
        need_replan = true;
    end

    if need_replan
        segment_active = true;
        segment_t0 = tk;
        segment_tau = tau_handover;

        d_plan = dHatRLS(:,k);

        % desired terminal state at handover time
        x_target_handover = target_state_from_roe(a_delta_alpha_target, segment_tau, w);

        % paper-style planning equation:
        % x_plan0 = Phi^{-1}(tau) [x_target - Gamma(tau) d_plan]
        Phi_seg = Phi_paper(segment_tau, w);
        Gamma_seg = Gamma_paper(segment_tau, w);

        % incremental correction version to avoid huge excursions
        x_free_end = Phi_seg*Xmeas(:,k) + Gamma_seg*d_plan;
        dx_needed = x_target_handover - x_free_end;

        dx_max = [20; 5; 20; 0.05; 0.02; 0.05];
        dx_cmd = min(max(dx_needed, -dx_max), dx_max);

        x_plan0 = Xmeas(:,k) + (Phi_seg \ dx_cmd);
    end

    % ------------------------------------------------------------
    % 10.6 Time-varying guidance reference
    % ------------------------------------------------------------
    s = tk - segment_t0;
    s = max(0, min(s, segment_tau));

    Phi_s = Phi_paper(s, w);
    Gamma_s = Gamma_paper(s, w);

    xGuid(:,k) = Phi_s*x_plan0 + Gamma_s*d_plan;

    % % ------------------------------------------------------------
    % % 10.7 PD tracking control
    % % ------------------------------------------------------------
    % pos_err = xGuid(1:3,k) - Xmeas(1:3,k);
    % vel_err = xGuid(4:6,k) - Xmeas(4:6,k);
    % 
    % pos_err_max = [20; 5; 20];
    % vel_err_max = [0.05; 0.02; 0.05];
    % 
    % pos_err = min(max(pos_err, -pos_err_max), pos_err_max);
    % vel_err = min(max(vel_err, -vel_err_max), vel_err_max);
    % 
    % u = Kp*pos_err + Kd*vel_err;
    % 
    % % component saturation
    % u = min(max(u, -umax), umax);
    % 
    % % norm saturation
    % if norm(u) > umax
    %     u = umax * u / norm(u);
    % end
    % 
    % uCmd(:,k) = u;
    % trackErr(k) = norm(xGuid(1:3,k) - X(1:3,k));

    % ------------------------------------------------------------
    % 10.8 Propagate truth
    % ------------------------------------------------------------
    X(:,k+1) = Ad*X(:,k) + Bd*(u + dTruth(:,k));
end

% finalize last-sample arrays
Xmeas(:,N) = X(:,N);
rhoHat(N) = rhoHat(N-1);
dHatRLS(:,N) = dHatRLS(:,N-1);
dHatFilt(:,N) = dHatFilt(:,N-1);
uCmd(:,N) = uCmd(:,N-1);
xGuid(:,N) = xGuid(:,N-1);
trackErr(N) = trackErr(N-1);

vrel_mag_truth = norm(X(4:6,N));
phiTruth(N) = 0.5*Cd*(Ac/mc)*(Vt + vrel_mag_truth)^2 ...
            - 0.5*Cd*(At/mt)*(Vt)^2;
dTruth(:,N) = [rhoTrue(N) * phiTruth(N); 0; 0];

%% 11) Plots - trajectory
figure();
plot3(X(2,:), X(3,:), X(1,:), 'LineWidth', 1.4);
grid on;
xlabel('y [m]');
ylabel('z [m]');
zlabel('x [m]');
title('Relative trajectory with guidance and control');
view(35,20);

figure();
plot(X(3,:), X(1,:), 'LineWidth', 1.5); hold on;
plot(xGuid(3,:), xGuid(1,:), '--', 'LineWidth', 1.2);
grid on;
xlabel('z [m]');
ylabel('x [m]');
title('x-z trajectory');
legend('Truth trajectory','Guidance reference','Location','best');

%% 12) Plots - drag estimation
figure();
plot(t/Torb, dTruth(1,:), 'k', 'LineWidth', 1.6); hold on;
plot(t/Torb, dHatFilt(1,:), 'r', 'LineWidth', 1.0);
plot(t/Torb, dHatRLS(1,:), 'b', 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('d_x [m/s^2]');
title('Differential drag estimation');
legend('Truth', 'DDOB filtered', 'DDOB + RLS', 'Location', 'best');

%% 13) Plots - density estimate
figure();
plot(t/Torb, rhoTrue, 'k', 'LineWidth', 1.6); hold on;
plot(t/Torb, rhoHat, 'b', 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('\rho [kg/m^3]');
title('Atmospheric density estimation (RLS)');
legend('Truth', '\rho-hat', 'Location', 'best');

%% 14) Plots - control and tracking
figure();
plot(t/Torb, vecnorm(uCmd,2,1), 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('|u| [m/s^2]');
title('Control acceleration magnitude');

figure();
plot(t/Torb, trackErr, 'LineWidth', 1.4);
grid on;
xlabel('Time [orbits]');
ylabel('Position tracking error [m]');
title('Guidance tracking error');

%% 15) Error summary
errRaw  = dHatRaw - dTruth;
errFilt = dHatFilt - dTruth;
errRLS  = dHatRLS - dTruth;

rmseRaw  = sqrt(mean(vecnorm(errRaw,2,1).^2));
rmseFilt = sqrt(mean(vecnorm(errFilt,2,1).^2));
rmseRLS  = sqrt(mean(vecnorm(errRLS,2,1).^2));
rhoRMSE  = sqrt(mean((rhoHat - rhoTrue).^2));
dV_used  = trapz(t, vecnorm(uCmd,2,1));

fprintf('\nDDOB raw RMSE      = %.8e m/s^2\n', rmseRaw);
fprintf('DDOB filt RMSE     = %.8e m/s^2\n', rmseFilt);
fprintf('DDOB + RLS RMSE    = %.8e m/s^2\n', rmseRLS);
fprintf('Density RMSE       = %.8e kg/m^3\n', rhoRMSE);
fprintf('Approx. control dV = %.8e m/s\n', dV_used);

%% ------------------------------------------------------------------------
% Local functions
% -------------------------------------------------------------------------

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

function Gamma = Gamma_paper(tau, w)
    c = cos(w*tau);
    s = sin(w*tau);

    Gamma = [ ...
        4*(1-c)/w^2 - 1.5*tau^2,  0,  2*(w*tau - s)/w^2;
        0,                        (1-c)/w^2, 0;
       -2*(w*tau - s)/w^2,        0,  (1-c)/w^2;
        4*s/w - 3*tau,            0,  2*(1-c)/w;
        0,                        s/w, 0;
       -2*(1-c)/w,                0,  s/w];
end

function x_target = target_state_from_roe(a_delta_alpha_target, tau, w)
    c = cos(w*tau);
    s = sin(w*tau);

    M35 = [ ...
         0,    1,      2*s,      -2*c,     0,       0;
         0,    0,      0,         0,      -s,       c;
        -1,    0,      c,         s,       0,       0;
       -1.5*w, 0,    2*w*c,     2*w*s,     0,       0;
         0,    0,      0,         0,    -w*c,    -w*s;
         0,    0,   -w*s,       w*c,       0,       0];

    x_target = M35 * a_delta_alpha_target;
end

function y = dt_placeholder(x)
    y = x;
end