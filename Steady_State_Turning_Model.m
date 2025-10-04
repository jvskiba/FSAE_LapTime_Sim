% turn_sim_fixed.m
% Adapted from user's original code. Keeps original structure where possible,
% fixes indexing/meshgrid/undefined functions and adds an iterative equilibrium
% solve for ax/ay. Also includes a simple linear tire model with saturation.

clear; clc; close all;

%% parameterize sweep
alpha_rl_vals = 0:0.1:12;   % sweep for rear-left slip angle (deg)
alpha_fl_vals = 0:0.1:12;   % sweep for front-left slip angle (deg)
[alpha_rl_grid, alpha_fl_grid] = meshgrid(alpha_rl_vals, alpha_fl_vals);

n_rl = numel(alpha_rl_vals);
n_fl = numel(alpha_fl_vals);

%% Pre-allocate
M_yaw = zeros(n_fl, n_rl);   % rows = alpha_fl index (i), cols = alpha_rl index (k)

%% Main nested loops (kept your loop style)
for k = 1:n_rl
    alpha_rl = alpha_rl_vals(k);
    for i = 1:n_fl
        alpha_fl = alpha_fl_vals(i);

        %% Inputs (kept from original)
        m = 220+68; % vehicle mass [kg]
        g = 9.81; % gravity [m/s^2]
        zcg = 0.300; % cg height [m]
        fWp = 47; % percent of weight on front axle [%]
        lWp = 50; % percent of weight on left side [%]
        P = 10; % tire pressure [psi] (passed to tire model but unused in stub)
        IA = 0; % tire inclination angle [deg]

        L = 1.530; % wheelbase [m]
        Tf = 1.368; % front track [m]
        Tr = 1.368; % rear track [m]
        Tav = (Tf+Tr)/2; % average track width

        kf = 373; % front roll stiffness [Nm/deg]
        kr = 392; % rear roll stiffness [Nm/deg]

        % steering inputs
        delta_fl = 20; % steering angle of FL [deg] (kept as in original)
        % original polynomial (kept as-is)
        delta_fr = 3.89 ... % steering angle of FR [deg]
            - 2.7e-18*delta_fl ...
            - 4.78e-05*delta_fl.^2 ...
            - 9.65e-23*delta_fl.^3 ...
            - 6.8e-10*delta_fl.^4 ...
            + 1.78e-26*delta_fl.^5 ...
            + 6.93e-14*delta_fl.^6 ...
            - 3.01e-30*delta_fl.^7 ...
            - 1.01e-17*delta_fl.^8 ...
            + 1.5e-34*delta_fl.^9 ...
            + 4.21e-22*delta_fl.^10;

        %% Intermediate geometry (kept)
        a = (100-fWp)/100 * L; % longitudinal distance from front axle to CG
        b = L - a;

        %% initial guesses for accelerations
        ax = 0;                % longitudinal accel [m/s^2]
        ay = 1*g;              % lateral accel guess (start at 1g to help convergence sometimes)

        %% initial corner weights (static)
        W_fl = m*g*b/L * lWp/100;
        W_fr = m*g*b/L * (100-lWp)/100;
        W_rl = m*g*a/L * lWp/100;
        W_rr = m*g*a/L * (100-lWp)/100;

        % error check (kept)
        W_tot = W_fl + W_fr + W_rl + W_rr;
        W_er = W_tot - m*g;
        if abs(W_er) >= .001*m*g
            fprintf('warning: initial corner weight mismatch (%.3f N)\n', W_er);
        end

        %% Iterative solution for steady-state forces/accelerations
        % We'll iterate until ax/ay converge (forces consistent with acceleration)
        max_iter = 200;
        tol_acc = 1e-4; % [m/s^2]
        prev_ax = ax; prev_ay = ay;

        for iter = 1:max_iter
            % ---- weight transfer based on current ax, ay ----
            d_fW = m * -ax * zcg / L;    % longitudinal transfer (N)
            d_lW = m * -ay * zcg / Tav;  % lateral transfer (N)

            % dynamic normal loads (kept formulas)
            N_fl = W_fl + (kf/(kf+kr)) * d_lW + d_fW;
            N_fr = W_fr - (kf/(kf+kr)) * d_lW + d_fW;
            N_rl = W_rl + (kr/(kf+kr)) * d_lW - d_fW;
            N_rr = W_rr - (kr/(kf+kr)) * d_lW - d_fW;

            % error check
            N_tot = N_fl + N_fr + N_rl + N_rr;
            N_er = N_tot - m*g;
            if abs(N_er) >= .01*m*g
                % not fatal, but report if big
                % fprintf('warning: dynamic normal loads mismatch (%.3f N)\n', N_er);
            end

            % ---- tire positions relative to CG (kept) ----
            p_fl = [a,  Tf/2];
            p_fr = [a, -Tf/2];
            p_rl = [-b, Tr/2];
            p_rr = [-b, -Tr/2];

            % ---- ICR intersection (kept method) ----
            % slopes of the front and rear tire force directions
            m_fl = tand(90 - (delta_fl - alpha_fl));
            m_rl = tand(90 + alpha_rl);

            b_fl = p_fl(2) - m_fl * p_fl(1);
            b_rl = p_rl(2) - m_rl * p_rl(1);

            % handle parallel/near-parallel lines (avoid division by nearly zero)
            denom = (m_fl - m_rl);
            if abs(denom) < 1e-6
                % fallback: set a large radius (near-straight)
                x_int = 1e6;
                y_int = m_fl * x_int + b_fl;
            else
                x_int = (b_rl - b_fl) / denom;
                y_int = m_fl * x_int + b_fl;
            end

            ICR = [x_int, y_int];
            R = sqrt(x_int^2 + y_int^2);

            % ---- inside slip angles (kept formulas) ----
            alpha_fr = (90 - atand((ICR(2) - p_fr(2)) / (ICR(1) - p_fr(1)))) - delta_fr;
            alpha_rr = 90 + atand((ICR(2) - p_rr(2)) / (ICR(1) - p_rr(1)));

            % ---- Tire forces via tire model ----
            % Note: UGAMS_TireModel is implemented below in this file.
            F_fl = UGAMS_TireModel(alpha_fl, N_fl, P, IA);
            F_fr = UGAMS_TireModel(alpha_fr, N_fr, P, IA);
            F_rl = UGAMS_TireModel(alpha_rl, N_rl, P, IA);
            F_rr = UGAMS_TireModel(alpha_rr, N_rr, P, IA);

            % ---- longitudinal forces due to steering geometry (kept) ----
            Fx_fl = -sind(delta_fl) * F_fl;
            Fx_fr = -sind(delta_fr) * F_fr;

            % rear longitudinal force allocation to satisfy steady-state sum
            % (kept from original)
            denom_Nrear = (N_rl + N_rr);
            if denom_Nrear == 0
                Fx_rl = 0;
            else
                Fx_rl = -(Fx_fl + Fx_fr) * N_rl / denom_Nrear;
            end
            Fx_rr = -(Fx_fl + Fx_fr) - Fx_rl;

            % ---- lateral forces (kept) ----
            Fy_fl = cosd(delta_fl) * F_fl;
            Fy_fr = cosd(delta_fr) * F_fr;
            Fy_rl = F_rl - Fx_rl;   % rear lateral uses remainder of vector
            Fy_rr = F_rr - Fx_rr;

            % ---- totals ----
            Fx_total = Fx_fl + Fx_fr + Fx_rl + Fx_rr;
            Fy_total = Fy_fl + Fy_fr + Fy_rl + Fy_rr;

            % ---- update accelerations ----
            ax_new = Fx_total / m;
            ay_new = Fy_total / m;

            % check convergence
            if (abs(ax_new - ax) < tol_acc) && (abs(ay_new - ay) < tol_acc)
                ax = ax_new;
                ay = ay_new;
                break;
            end

            % relaxation to improve stability
            relax = 0.5;
            ax = (1-relax)*ax + relax*ax_new;
            ay = (1-relax)*ay + relax*ay_new;
        end % iter

        % ---- Yaw moment (kept formulas) ----
        M_fl = -a*Fy_fl - Tf/2 * Fx_fl;
        M_fr = -a*Fy_fr + Tf/2 * Fx_fr;
        M_rl = b*Fy_rl - Tr/2 * Fx_rl;
        M_rr = b*Fy_rr + Tr/2 * Fx_rr;

        M_yaw(i, k) = M_fl + M_fr + M_rl + M_rr;

        % (optionally) store converged ax/ay if you want to analyze later
        % ax_store(i,k) = ax;
        % ay_store(i,k) = ay;
    end % alpha_fl
end % alpha_rl

%% Plotting results (kept and fixed grid variables)
tol = 1;   % Nm tolerance for near-zero yaw moment

figure('Name','Yaw Moment Surface','NumberTitle','off');
surf(alpha_rl_grid, alpha_fl_grid, M_yaw);
xlabel('\alpha_{rl} (deg)');
ylabel('\alpha_{fl} (deg)');
zlabel('M_{yaw} (Nm)');
title('Yaw Moment Surface');
shading interp;
colorbar;
hold on;

% Find points near zero yaw moment
[rows, cols] = find(abs(M_yaw) <= tol);

if ~isempty(rows)
    % Convert indices to actual coordinates
    alpha_rl_zero = alpha_rl_grid(sub2ind(size(M_yaw), rows, cols));
    alpha_fl_zero = alpha_fl_grid(sub2ind(size(M_yaw), rows, cols));
    M_yaw_zero    = M_yaw(sub2ind(size(M_yaw), rows, cols));

    % Overlay near-zero points in red
    plot3(alpha_rl_zero, alpha_fl_zero, M_yaw_zero, 'ro', ...
        'MarkerFaceColor','r', 'MarkerSize',4);
end

hold off;

%% Simple Tire Model function (linear + saturation)
% This is a stand-in for your UGAMS_TireModel. It accepts slip angle in degrees
% and returns a single scalar "magnitude" force (N) aligned with the tire
% resultant direction in the main script (so Fx/Fy decomposition is handled there).
%
% Replace this function with your real UGAMS_TireModel (Pacejka, lookup, etc.)
% when available.

function F = UGAMS_TireModel(alpha_deg, N, P, IA)
    % Parameters (tunable)
    mu = 1.3;             % friction coefficient (typical high-grip tire)
    C_alpha = 8000;       % cornering stiffness [N/rad] approximate
    % Convert slope to rad
    alpha = deg2rad(alpha_deg);

    % linear lateral force magnitude
    F_lin = C_alpha * alpha;  % N (signed)

    % saturation limit (friction circle approx - we only have lateral here)
    F_max = mu * N;

    % clip
    F = max(min(F_lin, F_max), -F_max);

    % If normal load is negative or non-physical, set zero
    if N <= 0
        F = 0;
    end
end
