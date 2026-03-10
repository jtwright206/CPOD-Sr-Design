%% yokota_2026_ddob_rls_rendezvous_demo.m
% Differential Drag Estimation (DDOB + RLS) + Autonomous Trajectory Correction
% Based on: Yokota et al., "Differential Drag Estimation and Autonomous Control
% for Non-Cooperative Rendezvous", AIAA SciTech 2026.
%
% Implements:
%  - CW/Hill relative dynamics (circular orbit) state-space (Eqs. 7–12)
%  - DDOB (Eq. 18) with a low-pass filter Q
%  - Differential drag model and RLS density estimator (Eqs. 20–26)
%  - Trajectory drift estimator x_target_center_hat (Eq. 27) + trigger (Eq. 28)
%  - Simple PD feedback with saturation (paper: PD + DOB + anti-windup saturation)
%
% Notes:
%  - This script is a compact "engineering implementation" of the paper's blocks,
%    not a perfect reproduction of their full simulator (J2, Harris-Priester, etc.).
%  - You can add J2 modeled disturbance a_d and/or a higher fidelity atmosphere later.

clear; clc; close all;

%% ---------------------- Parameters (Table 1-ish) -------------------------
% Gains/filters (Table 1 gives PD gains and Q cutoff; use as starting point)
Kp = 3e-2;          % [1/s^2] proportional gain (applied to position error)
Kd = 3e-1;          % [1/s]   derivative gain (applied to velocity error)
f_c_Q = 10;         % [Hz] low-pass cutoff for Q (Table 1)
lambda = 0.999;     % RLS forgetting factor (Table 1)
Tmax = 8;           % [N] thruster limit (Table 1), will convert to accel via m_c

% Spacecraft parameters (Table 1)
m_t = 3000;   A_t = 40;   Cd_t = 2.5;
m_c = 500;    A_c = 1;    Cd_c = 2.5;

% Orbit / CW parameters (choose an example circular orbit)
mu = 3.986004418e14;        % [m^3/s^2]
Re = 6378.137e3;            % [m]
alt = 550e3;                % [m] pick LEO altitude (paper uses LEO; exact alt not critical here)
r0  = Re + alt;
omega = sqrt(mu/r0^3);      % [rad/s]

% Simulation length: "4 revolutions" like their plots
Torb = 2*pi/omega;
Nrev = 4;
t_end = Nrev*Torb;

% Guidance / handover constraint (paper: x_target_center = -300 m, threshold 10 m)
x_target_center_nom = -300;     % [m]
dx_threshold = 10;              % [m]

% Handover time (choose the end-of-sim as "target time" or pick 1 rev ahead)
t_target = t_end;  % [s] nominal handover time

% Navigation noise (Table 1: 1%)
nav_noise_frac = 0.01;

%% ---------------------- Initial relative state ---------------------------
% Paper gives initial a*delta-alpha (Eq. 39). Here we just set a reasonable
% initial Hill state consistent with "behind & below" geometry.
x0  = 0; y0 = 0; z0 = 0;     % [m]
xd0 = 0;    yd0 = 0; zd0 = 0;      % [m/s]
x = [0; 0; 0; 0; 0; 0];

%% ---------------------- CW matrices (Eqs. 10–12) ------------------------
A = [ 0 0 0 1 0 0;
      0 0 0 0 1 0;
      0 0 0 0 0 1;
      0 0 0 0 0 2*omega;
      0 -omega^2 0 0 0 0;
      0 0 3*omega^2 -2*omega 0 0 ];

B = [ 0 0 0;
      0 0 0;
      0 0 0;
      1 0 0;
      0 1 0;
      0 0 1 ];

%% ---------------------- Discrete-time setup -----------------------------
dt = 0.5;                      % [s] integration timestep
t = (0:dt:t_end).';
Nt = numel(t);

% Discretize with simple forward Euler (fine for demo; replace with expm for higher fidelity)
Ad = eye(6) + A*dt;
Bd = B*dt;

%% ---------------------- Low-pass filter Q for DDOB/RLS -------------------
% Implement Q as 1st-order LPF per axis: y_dot = wc*(u - y)
wc = 2*pi*f_c_Q;
Q_state = zeros(3,1);   % filter state for DDOB output magnitude path (per-axis)
Q_phi   = 0;            % filter state for scalar regressor phi (for phase alignment)

% Helper inline LPF step
lpf1 = @(y,u) y + dt*wc*(u - y);

%% ---------------------- RLS state (Eqs. 21–26) ---------------------------
rho_hat = 1e-12;       % initial density guess [kg/m^3], rough LEO scale
P = 1e24;              % large initial covariance (scalar RLS)

%% ---------------------- Storage for plots -------------------------------
Xhist  = zeros(Nt,6);
Uhist  = zeros(Nt,3);
d_true = zeros(Nt,1);
d_ddob = zeros(Nt,1);
d_hat  = zeros(Nt,1);
rhoH   = zeros(Nt,1);
xcentH = zeros(Nt,1);
fb_on  = false(Nt,1);

%% ---------------------- "Truth" environment model ------------------------
% Use a simple time-varying density to mimic orbital variation (paper uses Harris-Priester).
% This gives the estimator something meaningful to track.
rho0 = 3e-12;                 % baseline density
rho_amp = 2e-12;              % variation
rho_true_fun = @(tt) rho0 + rho_amp*sin(2*pi*tt/Torb);

% Target airspeed magnitude (rough LEO circular speed)
Vt_mag = sqrt(mu/r0);         % ~7.6 km/s

% Differential drag acts primarily along-track (y-axis in Hill is often cross-track;
% but paper defines tangential drift in x; their oval is x-z. For simplicity, we inject
% differential drag along x (tangential-like in their x-z view). You can rotate as needed.
drag_axis = [1;0;0];          % apply along x for demo

%% ---------------------- Main simulation loop ----------------------------
for k = 1:Nt
    tt = t(k);

    % --- "Navigation" measurement with noise ---
    x_nav = x;
    x_nav(1:3) = x_nav(1:3) .* (1 + nav_noise_frac*randn(3,1));
    x_nav(4:6) = x_nav(4:6) .* (1 + nav_noise_frac*randn(3,1));

    % --- True differential drag acceleration (Eq. 20, simplified to magnitude) ---
    rho_true = rho_true_fun(tt);

    % Relative velocity magnitude approximation:
    v_rel = norm(x(4:6));  % [m/s] relative speed
    Vc_mag = Vt_mag + v_rel;  % crude approximation (paper uses |Vt + v|)

    coeff_c = 0.5*Cd_c*A_c/m_c;
    coeff_t = 0.5*Cd_t*A_t/m_t;

    dmag_true = rho_true*( coeff_c*Vc_mag^2 - coeff_t*Vt_mag^2 ); % [m/s^2] magnitude
    dvec_true = dmag_true * drag_axis;  % [m/s^2]
    d_true(k) = dmag_true;

    % --- Controller (PD) computes commanded accel, but FB is enabled/disabled ---
    % Guidance reference: keep nominal "oval center" (very simplified). We just
    % hold a desired position in x-z around the nominal center.
    x_ref = [x_target_center_nom; 0; 0];  % [m] desired relative position (toy)
    v_ref = [0;0;0];

    pos_err = x_ref - x_nav(1:3);
    vel_err = v_ref - x_nav(4:6);

    a_cmd_fb = Kp*pos_err + Kd*vel_err;  % [m/s^2]

    % Convert accel command to thrust-limited accel:
    a_max = Tmax/m_c;  % [m/s^2]
    a_cmd_fb = max(min(a_cmd_fb, a_max), -a_max);

    % --- DDOB (Eq. 18 specialized to circular orbit: k = sqrt(omega), omega_dot=0) ---
    % Need x_ddot estimate from nav: use finite difference on velocity (noisy!)
    if k == 1
        v_prev = x_nav(4:6);
        a_nav_ddot = zeros(3,1);
    else
        a_nav_ddot = (x_nav(4:6) - v_prev)/dt;
        v_prev = x_nav(4:6);
    end

    % Predicted dynamics term (from Eq. 5 with omega_dot=0, k=omega^(1/2) for circular):
    % Using simplified Eq. 7 form directly:
    x_pred = x_nav(1); y_pred = x_nav(2); z_pred = x_nav(3);
    xd_pred = x_nav(4); yd_pred = x_nav(5); zd_pred = x_nav(6);

    f_dyn = [ 2*omega*zd_pred;
              -omega^2*y_pred;
              -2*omega*xd_pred + 3*omega^2*z_pred ];

    a_d_model = [0;0;0]; % modeled disturbance (e.g., J2) not included in this demo

    % DDOB raw estimate (before Q):
    d_ddob_raw = (a_nav_ddot - a_cmd_fb) - f_dyn - a_d_model;  % ~ d (plus noise)

    % Apply Q (LPF) per-axis:
    for i = 1:3
        Q_state(i) = lpf1(Q_state(i), d_ddob_raw(i));
    end
    dvec_ddob = Q_state;                % filtered DDOB vector estimate
    dmag_ddob = norm(dvec_ddob);        % use magnitude like Eq. 22
    d_ddob(k) = dmag_ddob;

    % --- RLS density estimator (Eqs. 21–26) ---
    % Regressor (Eq. 23): phi = Q * | 0.5*Cd_c*A_c/m_c*(Vt+v)^2 - 0.5*Cd_t*A_t/m_t*Vt^2 |
    phi_raw = abs( coeff_c*(Vt_mag + norm(x_nav(4:6)))^2 - coeff_t*Vt_mag^2 );
    Q_phi = lpf1(Q_phi, phi_raw);   % include Q for phase alignment
    phi = Q_phi;

    y = dmag_ddob;                  % Eq. 22

    % Scalar RLS update (stable form)
    % K = P*phi / (lambda + phi^2*P)
    denom = lambda + (phi^2)*P;
    K_rls = (P*phi)/denom;
    eps   = y - phi*rho_hat;

    rho_hat = rho_hat + K_rls*eps;
    P = (1/lambda)*(P - K_rls*phi*P);

    % Estimated differential drag (Eq. 26)
    dmag_hat = rho_hat * phi_raw;    % (use unfiltered phi_raw for fast v update)
    d_hat(k) = dmag_hat;
    rhoH(k)  = rho_hat;

    % --- Trajectory drift estimator for x_target_center_hat (Eq. 27) ---
    % Use x,z and their rates; and ax_hat, az_hat = (u + mean(d_hat) over a window).
    % Here: ax_hat approx = a_cmd + d_hat along x, az_hat = a_cmd_z (no drag axis on z).
    tau = (t_target - tt);
    if tau < 0, tau = 0; end

    % Form ax_hat, az_hat:
    ax_hat = a_cmd_fb(1) + dmag_hat*drag_axis(1);
    az_hat = a_cmd_fb(3) + dmag_hat*drag_axis(3);

    % Eq. (27) coefficients:
    % xcenter_hat = [1, 6*omega*tau, -3*tau^2*omega] * [x; z; xdot; zdot] + [4/omega^2 - 3/2*tau^2, 2/omega*tau]*[ax_hat; az_hat]
    xcent_hat = [1, 6*omega*tau, -3*(tau^2)*omega, 2/omega] * [x_nav(1); x_nav(3); x_nav(4); x_nav(6)] ...
              + [4/omega^2 - 1.5*(tau^2), 2/omega*tau] * [ax_hat; az_hat];

    xcentH(k) = xcent_hat;

    % Feedback ON/OFF trigger (Eq. 28)
    fb_on(k) = (abs(x_target_center_nom - xcent_hat) > dx_threshold);

    % Apply FB gating: if not enabled, command is zero (coast) in this demo
    if fb_on(k)
        a_cmd = a_cmd_fb;
    else
        a_cmd = zeros(3,1);
    end

    % Store control
    Uhist(k,:) = a_cmd(:).';

    % --- Propagate true system ---
    % x_{k+1} = Ad*x + Bd*(u + d_true)
    u_total = a_cmd + dvec_true;   % control accel + true drag accel
    x = Ad*x + Bd*u_total;

    Xhist(k,:) = x.';
end

%% ---------------------- Plots (similar to paper figures) -----------------
figure; plot3(Xhist(:,1), Xhist(:,3), Xhist(:,2)); grid on;
xlabel('x [m]'); ylabel('z [m]'); zlabel('y [m]');
title('Relative trajectory (Hill)');

figure; plot(t/Torb, d_true, 'k', t/Torb, d_ddob, 'r', t/Torb, d_hat, 'b'); grid on;
xlabel('Time [rev]'); ylabel('Differential drag magnitude [m/s^2]');
legend('True', 'DDOB (|Q d|)', 'Proposed (DDOB+RLS)', 'Location','best');
title('Differential drag estimation (paper Fig. 8 concept)');

figure; plot(t/Torb, rho_true_fun(t), 'k', t/Torb, rhoH, 'b'); grid on;
xlabel('Time [rev]'); ylabel('\rho [kg/m^3]');
legend('True', '\rho-hat (RLS)', 'Location','best');
title('Atmospheric density estimation (RLS)');

figure; plot(t/Torb, xcentH, 'b', 'LineWidth',1.2); hold on; grid on;
yline(x_target_center_nom,'k--','Nominal');
yline(x_target_center_nom+dx_threshold,'r--','Threshold');
yline(x_target_center_nom-dx_threshold,'r--','Threshold');
xlabel('Time [rev]'); ylabel('Estimated x_{center}(t_{target}) [m]');
title('Trajectory drift estimator / FB trigger (paper Fig. 9 concept)');

figure; stairs(t/Torb, fb_on, 'LineWidth',1.2); grid on;
xlabel('Time [rev]'); ylabel('FB ON (1=yes)');
title('Autonomous feedback enable logic');

figure; plot(t/Torb, vecnorm(Uhist,2,2)); grid on;
xlabel('Time [rev]'); ylabel('|a_{cmd}| [m/s^2]');
title('Control effort (accel command magnitude)');

%% ---------------------- End of script ------------------------------------