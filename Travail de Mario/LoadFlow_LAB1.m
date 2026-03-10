%% Load Flow Solver
%
% Schéma :
%
%    GCP
%    PM5        PM4   PM3   PM2   PM1
%     |          |     |     |     |
%     0---- ---- 0---- 0---- 0---- 0
%       Z5   Z4    Z3    Z2    Z1
%
% Node ordering (slack bus = node 1):
%   Node 1 = GCP (PM5) - Slack bus
%   Node 2 = PM4
%   Node 3 = PM3
%   Node 4 = PM2
%   Node 5 = PM1

clear; clc; close all;

%% =============================================
%  Base values
%  =============================================
Vbase = 230;          % [V] - nominal phase voltage
Sbase = 1e5;          % [VA] = 100 kVA
Zbase = Vbase^2 / Sbase;  % = 0.529 Ohm

%% =============================================
%  Line parameters (in PU)
%  =============================================
f = 50; % Hz
omega = 2 * pi * f;

% Z5 (buried cable, 60 m) - from datasheet:
%   resistance: 3.08 Ohm/km
%   reactance: 0.09 Ohm/km
%   Length = 60 m
L_Z5 = 0.060; % km
R_Z5_SI = 3.08 * L_Z5;   % Ohm
X_Z5_SI = 0.09 * L_Z5;   % Ohm
Z5_SI = R_Z5_SI + 1j * X_Z5_SI;

% Z4 to Z1 (identical lines)
%   R = 160 mOhm, L = 460 µH
R_line_SI = 160e-3;              % Ohm
X_line_SI = omega * 460e-6;      % Ohm
Z_line_SI = R_line_SI + 1j * X_line_SI;

% Convert to per unit
Z5_pu = Z5_SI / Zbase;
Z4_pu = Z_line_SI / Zbase;
Z3_pu = Z_line_SI / Zbase;
Z2_pu = Z_line_SI / Zbase;
Z1_pu = Z_line_SI / Zbase;

%% =============================================
%  Bus Admittance Matrix Y (5x5)
%  =============================================
% Topology (radial):
%   Branch 1: Node 1 -- Node 2 (Z5)
%   Branch 2: Node 2 -- Node 3 (Z4)
%   Branch 3: Node 3 -- Node 4 (Z3)
%   Branch 4: Node 4 -- Node 5 (Z2)
%   (Note: Z1 in the topology diagram connects PM2-PM1 = Node 4-Node 5,
%    Z2 connects PM3-PM2 = Node 3-Node 4, etc.
%    But since Z1=Z2=Z3=Z4 are identical, the mapping doesn't matter)

n_nodes = 5;

% Line admittances
y5 = 1 / Z5_pu;  % Branch: node 1 - node 2
y4 = 1 / Z4_pu;  % Branch: node 2 - node 3
y3 = 1 / Z3_pu;  % Branch: node 3 - node 4
y2 = 1 / Z2_pu;  % Branch: node 4 - node 5

% Build Y bus matrix using standard rules:
%   Diagonal Yii = sum of admittances connected to node i
%   Off-diagonal Yij = -admittance between nodes i and j
Y = zeros(n_nodes, n_nodes);

% Off-diagonal terms
Y(1,2) = -y5;  Y(2,1) = -y5;
Y(2,3) = -y4;  Y(3,2) = -y4;
Y(3,4) = -y3;  Y(4,3) = -y3;
Y(4,5) = -y2;  Y(5,4) = -y2;

% Diagonal terms (no shunt elements)
Y(1,1) = y5;
Y(2,2) = y5 + y4;
Y(3,3) = y4 + y3;
Y(4,4) = y3 + y2;
Y(5,5) = y2

%% =============================================
%  Define nodal power injections (in pu)
%  =============================================
% Convention: positive = generation, negative = consumption
%
% Power injections in Watts and VAR, then convert to pu
% Negative P = consumption, Positive P = generation
P2 = 0;          % PM4: 
Q2 = 0;          % PM4:
P3 = -5000;      % PM3: (5 kW load)
Q3 = -2000;      % PM3: (2 kVAR inductive)
P4 = -3000;      % PM2: 
Q4 = -1000;      % PM2: 
P5 = -4000;      % PM1: 
Q5 = -1500;      % PM1: 

% Convert to per unit
S_inj = ([P2+1j*Q2; P3+1j*Q3; P4+1j*Q4; P5+1j*Q5]) / Sbase

% Vectors of P and Q for non-slack nodes (nodes 2 to 5)
P_spec = real(S_inj);  % (4x1)
Q_spec = imag(S_inj);  % (4x1)


%% =============================================
%  Newton-Raphson Load Flow Solver
%  =============================================

% Initial guess: flat start (all voltages = 1 pu, angle = 0)
Vslack = 1 * exp(1j * 0);

V = [Vslack; Vslack; Vslack; Vslack; Vslack];  % 5x1 complex

% Solver parameters
tol = 1e-6;
max_iter = 100;
converged = false;

% Function to compute complex power injections from voltages
S_calc = @(V) diag(V) * conj(Y * V);

for iter = 1:max_iter
    
    % Compute power injections at current voltage estimate
    S_computed = S_calc(V);
    
    % Extract computed P and Q for non-slack buses (nodes 2:5)
    P_comp = real(S_computed(2:end));
    Q_comp = imag(S_computed(2:end));
    
    % Power mismatch vector: b = [P_spec - P_comp; Q_spec - Q_comp]
    mismatch = [P_spec - P_comp; Q_spec - Q_comp];
    
    % Check convergence
    max_mismatch = max(abs(mismatch));    
    if max_mismatch < tol
        converged = true;
        fprintf('Converged after %d iterations!\n', iter);
        break;
    end
    
    % Extract voltage magnitudes and angles for all nodes
    V_abs = abs(V);
    V_angle = angle(V);
    
    % Compute Jacobian (using the provided function)
    J = make_jacobian_matrix(V_abs, V_angle, Y, n_nodes);
    
    % Newton update: delta = J \ mismatch
    delta = J \ mismatch;
    
    % delta is structured as [delta_v; delta_angle] for non-slack nodes
    n_pq = n_nodes - 1;  % number of non-slack nodes = 4
    delta_v     = delta(1:n_pq);
    delta_angle = delta(n_pq+1:end);
    
    % Update non-slack bus voltages
    V_abs_new   = V_abs(2:end)   + delta_v;
    V_angle_new = V_angle(2:end) + delta_angle;
    
    % Reconstruct complex voltages for non-slack buses
    V_no_slack_new = V_abs_new .* exp(1j * V_angle_new);
    
    % Reconstruct full voltage vector
    V = [Vslack; V_no_slack_new];
end

if ~converged
    fprintf('WARNING: Did not converge within %d iterations!\n', max_iter);
end

%% =============================================
%  Results
%  =============================================
fprintf('\n=== Load Flow Results ===\n');
fprintf('%-8s %-14s %-14s %-14s\n', 'Node', '|V| (pu)', '|V| (V)', 'Angle (deg)');
fprintf('%-8s %-14s %-14s %-14s\n', '----', '--------', '--------', '-----------');
for k = 1:n_nodes
    fprintf('%-8d %-14.6f %-14.2f %-14.4f\n', k, abs(V(k)), abs(V(k))*Vbase, rad2deg(angle(V(k))));
end

% Compute slack bus power injection
S_final = S_calc(V);
P_slack = real(S_final(1)) * Sbase;
Q_slack = imag(S_final(1)) * Sbase;
fprintf('\nSlack bus power: P = %.2f W, Q = %.2f VAR\n', P_slack, Q_slack);

%% =============================================
%  Plot results
%  =============================================
figure('Name', 'Load Flow Results', 'Position', [100 100 900 400]);

subplot(1,2,1);
bar(1:n_nodes, abs(V));
xlabel('Node');
ylabel('|V| (pu)');
title('Voltage Magnitude');
ylim([min(abs(V))*0.995, max(abs(V))*1.005]);
grid on;
set(gca, 'XTickLabel', {'GCP','PM4','PM3','PM2','PM1'});

subplot(1,2,2);
bar(1:n_nodes, rad2deg(angle(V)));
xlabel('Node');
ylabel('Angle (degrees)');
title('Voltage Angle');
grid on;
set(gca, 'XTickLabel', {'GCP','PM4','PM3','PM2','PM1'});

%% =============================================
%  Line currents
%  =============================================
fprintf('\n=== Line Currents ===\n');
fprintf('%-12s %-20s %-14s\n', 'Line', 'Current (pu)', '|I| (pu)');
fprintf('%-12s %-20s %-14s\n', '----', '------------', '--------');

% Primitive admittance matrix (diagonal with line admittances)
% Current in branch from node i to node j: I_ij = y_ij * (Vi - Vj)
branches = {1, 2, y5, 'Z5 (1->2)';
            2, 3, y4, 'Z4 (2->3)';
            3, 4, y3, 'Z3 (3->4)';
            4, 5, y2, 'Z2 (4->5)'};

for b = 1:size(branches,1)
    from = branches{b,1};
    to   = branches{b,2};
    y_br = branches{b,3};
    name = branches{b,4};
    I_branch = y_br * (V(from) - V(to));
    fprintf('%-12s %.6f + j%.6f   |I| = %.6f pu\n', ...
        name, real(I_branch), imag(I_branch), abs(I_branch));
end