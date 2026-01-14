%% Helipcopter Blade Vibration Analysis
% Author: Christian Gonzalez Perez
% Date: 2026-01

clear; clc; close all;

%% 1. Parameters Definition
chi = 30; % ID card parameter
R = 6.25; % Blade Radius [m]
m0 = 6.05e-2; % Linear mass [kg/m]
E = 5e10; % Young Modulus [N/m^2]
RPM = 250 + (chi / 10);
Omega_val = RPM * (2 * pi / 60); % Rotation speed [rad/s]

sigma = 0.05; % Solidity
Nb = 4;       % Assuming 4 blades
c = (sigma * pi * R) / Nb; % Blade chord [m]
t = 0.12 * c; % NACA 0012 profile thickness
I0 = 0.04 * c * t^3; % Area Moment of Inertia for aerodynamic section
EI0 = E * I0; % Structural stiffness [N*m^2]

syms x y;
N_modes = 4; % Number of assumed modes
x_grid = linspace(0, 1, 101);   % Discretization of the blade

%% 2. Cases Setup
distributions= {'Uniform', 'Non_Uniform'};
bcTypes = {'hinged', 'cantilevered'};
rotationSpeeds = [0, Omega_val]; % 0 = Non-rotating, Omega_val = Rotating

results = struct();

%% 3. Execution Loop
for d = 1:length(distributions)
    dist = distributions{d};
    fprintf('Case: %s\n', dist);

    % Define symbolic distributions
    if strcmp(dist, 'Uniform')
        m_dist = m0;
        EI_dist = EI0;
    else
        m_dist = 4 * m0 * (1 - x);
        EI_dist = 4 * EI0 * (1 - x);
    end

    for b = 1:length(bcTypes)
        bc = bcTypes{b};
        gamma = buildGamma(N_modes, bc);

        for rot = 1:length(rotationSpeeds)
            Omega = rotationSpeeds(rot);
            rotLabel = 'NonRotating';
            if Omega > 0
                rotLabel = 'Rotating';
            end

            fprintf('BC: %s, Rotation speed: %d rad/s \n', bc, Omega);

            % Build Centrifugal Tension G(x)
            G = buildG(m_dist, Omega, R);

            % Assemble M and K
            [M_sym, K_sym] = assembleMK(gamma, m_dist, EI_dist, G, R);

            % Convert to numeric for eigenvalue solving
            M_num = double(M_sym);
            K_num = double(K_sym);

            % Solve generalized eigenvalue problem
            [vecs, vals] = eig(K_num, M_num); % Eigenvectors in columns
            omega2 = diag(vals); % Natural frequencies squared

            [omega2_sorted, idx] = sort(real(omega2)); % Sort frequencies
            vecs_sorted = vecs(:, idx); % Sort the modes

            selected_omega = sqrt(omega2_sorted(1:N_modes));
            selected_vecs = vecs_sorted(:, 1:N_modes);

            % Store results
            caseID = sprintf('%s_%s_%s', dist, bc, rotLabel);
            results.(caseID).frequencies = selected_omega;
            results.(caseID).lambda = selected_omega / (Omega + eps); % Avoid div by zero

            % Reconstruct Mode Shapes
            modes = zeros(N_modes, length(x_grid));
            for i = 1:N_modes
                v = selected_vecs(:, i);
                mode_sym = 0;
                for j = 1:N_modes
                    mode_sym = mode_sym + v(j) * gamma(j);
                end
                modes(i, :) = double(subs(mode_sym, x, x_grid));
                % Normalize: max absolute value = 1 and force positive
                modes(i, :) = modes(i, :) / max(abs(modes(i, :)));
                if modes(i, end) < 0
                    modes(i, :) = -modes(i, :);
                end
            end
            results.(caseID).modes = modes;
        end
    end
end

%% 4. Visualization
% Figure 1: Uniform Case - Hinged
plotComparative(results, 'Uniform', 'hinged', x_grid);
% Figure 2: Uniform Case - Cantilevered
plotComparative(results, 'Uniform', 'cantilevered', x_grid);
% Figure 3: Non-Uniform Case - Hinged
plotComparative(results, 'Non_Uniform', 'hinged', x_grid);
% Figure 4: Non-Uniform Case - Cantilevered
plotComparative(results, 'Non_Uniform', 'cantilevered', x_grid);
fprintf('\nAnalysis Complete.\n');

%% Helper: Plotting Function
function plotComparative(results, dist, bc, x_grid)
caseRot = sprintf('%s_%s_Rotating', dist, bc);
caseNon = sprintf('%s_%s_NonRotating', dist, bc);

figure('Name', sprintf('%s - %s', dist, bc), 'NumberTitle', 'off');
N_modes = 4; % Number of modes to plot

for i = 1:N_modes
    subplot(2, 2, i);
    plot(x_grid, results.(caseNon).modes(i,:), 'b--', 'LineWidth', 2);
    hold on;
    plot(x_grid, results.(caseRot).modes(i,:), 'r-', 'LineWidth', 1);
    grid on;
    title(sprintf('Mode %d', i));
    legend(sprintf('f_0=%.2f Hz', results.(caseNon).frequencies(i)/(2*pi)), ...
        sprintf('f_{rot}=%.2f Hz', results.(caseRot).frequencies(i)/(2*pi)), ...
        'Location', 'best');
    if i > 2, xlabel('x/R'); end
    if mod(i, 2) ~= 0, ylabel('Deflection'); end
end
sgtitle(sprintf('%s Blade - %s', dist, bc),'Interpreter','none');
end
%% Centrifugal Tension Function
function G = buildG(m_sym, Omega, R)
% buildG: Calculates centrifugal tension G(x)
% m_sym: symbolic mass distribution m(x)
% Omega, R: scalars

syms x y;
G = Omega^2 * R^2 * int(subs(m_sym, x, y) * y, y, x, 1);
end