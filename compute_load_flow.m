%% Structure of the grid lab
%
% Schéma :
%
%   GCP
%   PM5         PM4    PM3    PM2    PM1
%    |     a     |  b   |   c  |   d  |
%    1 ---- ---- 2 ---- 3 ---- 4 ---- 5
%       Z5   Z4     Z3     Z2     Z1
%
% Node ordering (slack bus = node 1):
%   Node 1 = GCP (PM5) - Slack bus
%   Node 2 = PM4
%   Node 3 = PM3
%   Node 4 = PM2
%   Node 5 = PM1

clear; clc; close all;

%% Definition of the base values
Vbase = 230; % 230 V (nominal phase voltage)
Sbase = 1e5; % 100 kVA
Zbase = Vbase^2/Sbase; % = 0.529 Ohm

%% Computation of the line parameters in pu

% Computation of line paramaters in Ohm
Z_5_SI = (3.08 + 0.09i) * 60 / 1000; % length x (linear resistance and linear reactance)
Z_others_SI = 0.16 + 460 * 1e-6i * 2 * pi * 50;

% Conversion to pu units
Z_5_pu = Z_5_SI / Zbase;
Z_others_pu = Z_others_SI / Zbase;

% Assigmnent of the value to the other impedances
Z_1_pu = Z_others_pu;
Z_2_pu = Z_others_pu;
Z_3_pu = Z_others_pu;
Z_4_pu = Z_others_pu;

fprintf("\nThe impedance of the GCP is %f %+fi\n", real(Z_5_pu), imag(Z_5_pu));
fprintf("\nThe impedance of the others lines are %f %+fi\n", real(Z_others_pu), imag(Z_others_pu));

% Computation of the admittances
Y_a = 1/(Z_5_pu + Z_4_pu);
Y_b = 1/Z_3_pu;
Y_c = 1/Z_2_pu;
Y_d = 1/Z_1_pu;

%% Definition of the branch to node matrix
A = [ 1, -1,  0,  0,  0;
      0,  1, -1,  0,  0;
      0,  0,  1, -1,  0;
      0,  0,  0,  1, -1];

% Computation of the number of nodes and of lines
n_lines = size(A, 1); % Extract the 1st dimension
n_nodes = size(A, 2); % Extract the 2nd dimension

fprintf("\nThe number of lines is %i\n", n_lines);
fprintf("\nThe number of nodes is %i\n", n_nodes);

%% Definition of the primitive and bus admittance matrix

% Definition of the primitive admittance matrix
Y_PR = [Y_a,   0,   0,   0;
          0, Y_b,   0,   0;
          0,   0, Y_c,   0;
          0,   0,   0, Y_d];

% Definition of the bus admittance matrix (with and without the slackbus)
Y_bus = transpose(A) * Y_PR * A;

Y_bus_no_slackbus = Y_bus(2:end, 2:end);

fprintf("\nThe bus admittance matrix is \n");
disp(Y_bus);

%% Initialization of the vector of nodal power injections (excludes the slackbus)

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
S_spec = ([P2+1j*Q2; P3+1j*Q3; P4+1j*Q4; P5+1j*Q5]) / Sbase;

fprintf("\nThe vector of nodal power injections in pu without the slackbus contains\n");
disp(S_spec);

%% Initialization of the voltages

% Initial guess for the voltage at the slackbus: flat start (all voltages = 1 pu, angle = 0)
V_slack = 1 * exp(1i*0);

% Initial guess for the vector of nodal voltages
V = V_slack * ones(n_nodes, 1);

fprintf("\nThe vector of nodal voltages in pu contains\n");
disp(V);

%% Application of the Newton method

% Definition of the stop criteria
tol = 1e-6;
max_iter = 300;
converged = false;

% Definition of voltage and power helper functions
to_complex = @(x) (x(1:numel(x)/2) .* exp(1i.*x(numel(x)/2+1:end)));
to_polar = @(z) [abs(z(:)); angle(z(:))];

calc_S = @(V, Y_bus) diag(V) * conj(Y_bus * V);
stack_PQ = @(S) vertcat(real(S), imag(S));

% Exclusion of the slackbus node
V_pq = V(2:end);
x = to_polar(V_pq);

% Use full V vector and full Y_bus to match the Jacobian computation
S_calc = calc_S(V, Y_bus);
S_calc_pq = S_calc(2:end);

mismatch = S_spec - S_calc_pq;

mismatch_pq = stack_PQ(mismatch);

% Computation of the error
error = norm(mismatch);

% Initialization
iter = 0;

while iter < max_iter && error > tol

    % Increment the number of iterations
    iter = iter + 1;

    % Calculate Jacobian (use the provided helper file)
    J = make_jacobian_matrix(abs(V), angle(V), Y_bus, n_nodes);

    % Perform Newton update step
    dx = J\(mismatch_pq); % Solve J * delta_x = mismatch

    % Update state vector
    x = x + dx;

    % Reconstruct PQ bus voltages
    V_pq = to_complex(x);

    % Reconstruct the entire voltage vector
    V = vertcat(V_slack, V_pq);

    % Compute the new approximated apparent power in complex form
    % Use full V vector and full Y_bus to match the Jacobian computation
    S_calc = calc_S(V, Y_bus);
    S_calc_pq = S_calc(2:end);

    % Computation of the error
    mismatch = S_spec - S_calc_pq;
    mismatch_pq = stack_PQ(mismatch);

    error = norm(mismatch);

    % Print some intermediate results
    if mod(iter, 10) == 1

        fprintf("------ ITERATION n°%i -------\n", iter);
        
        fprintf("The error is %f for the following power difference\n\n", error);
        disp(mismatch);

        fprintf("And the following voltage value \n\n");
        disp(V_pq);
    end

end

%% Display final results
fprintf('\n=== Load Flow Results ===\n');
fprintf('Converged after %d iterations!\n', iter);
fprintf('Final error: %.2e\n\n', error);

fprintf('%-8s %-14s %-14s %-14s\n', 'Node', '|V| (pu)', '|V| (V)', 'Angle (deg)');
fprintf('%-8s %-14s %-14s %-14s\n', '----', '--------', '--------', '-----------');
fprintf('%-8d %-14.6f %-14.2f %-14.4f\n', 1, abs(V_slack), abs(V_slack)*Vbase, rad2deg(angle(V_slack)));
for k = 2:n_nodes
    fprintf('%-8d %-14.6f %-14.2f %-14.4f\n', k, abs(V(k)), abs(V(k))*Vbase, rad2deg(angle(V(k))));
end

% Compute slack bus power injection
S_final = calc_S(V, Y_bus);
P_slack = real(S_final(1)) * Sbase;
Q_slack = imag(S_final(1)) * Sbase;
fprintf('\nSlack bus power: P = %.2f W, Q = %.2f VAR\n', P_slack, Q_slack);

%% =============================================
%  Plot results
%  =============================================
figure('Name', 'Load Flow Results', 'Position', [100 100 1200 400]);

% Subplot 1: Voltage Magnitude and Angle
subplot(1,2,1);
yyaxis left;
bar(1:n_nodes, abs(V)*Vbase);
ylabel('Magnitude [V]');
ylim([min(abs(V))*0.995, max(abs(V))*1.005]*Vbase);

yyaxis right;
plot(1:n_nodes, rad2deg(angle(V)), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0.85 0.33 0.1]);
ylabel('Angle [deg]');
grid on;
title('Voltage Profile');
xlabel('Node');
set(gca, 'XTickLabel', {'GCP','PM4','PM3','PM2','PM1'});
legend({'Voltage Magnitude', 'Voltage Angle'}, 'Location', 'best');

% Subplot 2: Power Injections
subplot(1,2,2);
% Reconstruct power injections from calculated values
P_nodes = real(S_final) * Sbase;  % Convert to Watts
Q_nodes = imag(S_final) * Sbase;  % Convert to VAR

% Create grouped bar chart
b = bar([P_nodes, Q_nodes]);
b(1).FaceColor = [0.2 0.6 0.8];  % Active power - blue
b(2).FaceColor = [0.8 0.4 0.2];  % Reactive power - orange

xlabel('Node');
ylabel('Power [W, VAR]');
title('Power Injections');
grid on;
set(gca, 'XTick', 1:n_nodes);
set(gca, 'XTickLabel', {'GCP','PM4','PM3','PM2','PM1'});
legend({'Active Power (P)', 'Reactive Power (Q)'}, 'Location', 'best');
ylim([min([0; P_nodes; Q_nodes]) * 1.05, max([0; P_nodes; Q_nodes]) * 1.05]);
