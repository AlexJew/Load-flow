%% ========================================================================
%  LOAD FLOW par la méthode de Newton-Raphson
%  CPRe - Load Flow
%  ========================================================================
%  Ce script implémente un load flow sur le réseau 4 bus du cours.
%
%  RÉSEAU (cf. schéma du cours) :
%
%              ②
%             / | \
%        Yb  / Yd  \ Ye
%           /   |   \
%    ③ ---Yc--- ① ---Yf--- ④
%     |                      |
%    Ya                     Yg
%     |                      |
%    ⓪ (référence/terre)   ⓪
%
%  Branches :
%    Ya : Bus 3 -> Bus 0 (shunt terre, impedance interne generateur 3)
%    Yb : Bus 3 -> Bus 2 (ligne)
%    Yc : Bus 3 -> Bus 1 (ligne)
%    Yd : Bus 2 -> Bus 1 (ligne)
%    Ye : Bus 2 -> Bus 4 (ligne)
%    Yf : Bus 1 -> Bus 4 (ligne)
%    Yg : Bus 4 -> Bus 0 (shunt terre, impedance interne generateur 4)
%
%  Bus 1 : Slack bus (V = 1.0 pu, angle = 0)
%  Bus 2, 3, 4 : Bus PQ
%
%  METHODE DE NEWTON-RAPHSON :
%    1. Initialiser x = [v; phi] (no load : v=1 pu, phi=0)
%    2. Calculer g(x) = mismatch de puissance
%    3. Calculer le Jacobien J
%    4. Resoudre J * delta = -g(x)
%    5. x_new = x_old + delta
%    6. Repeter jusqu'a ||g(x)||_inf <= tolerance
%  ========================================================================

clear; clc; close all;

%% ========================================================================
%  ETAPE 1 : DEFINITION DU RESEAU
%  ========================================================================

N = 4;          % Nombre de bus (1, 2, 3, 4)
slack_bus = 1;  % Bus 1 est le slack bus

% --- Impedances des branches (en pu) ---
Za = 0.03 + 1i*0.20;   % Bus 3 -> terre (impedance interne generateur 3)
Zb = 0.01 + 1i*0.10;   % Bus 3 -> Bus 2 (ligne)
Zc = 0.02 + 1i*0.15;   % Bus 3 -> Bus 1 (ligne)
Zd = 0.015 + 1i*0.12;  % Bus 2 -> Bus 1 (ligne)
Ze = 0.01 + 1i*0.08;   % Bus 2 -> Bus 4 (ligne)
Zf = 0.02 + 1i*0.10;   % Bus 1 -> Bus 4 (ligne)
Zg = 0.025 + 1i*0.18;  % Bus 4 -> terre (impedance interne generateur 4)

% --- Injections de puissance connues (en pu) ---
% Bus 1 (slack) : P et Q non utilises
% Bus 2 : charge de 1.0 + j0.4 pu
% Bus 3 : generation de 0.8 + j0.3 pu
% Bus 4 : generation de 0.5 + j0.2 pu
P_spec = [0;  -1.0;  +0.8;  +0.5];
Q_spec = [0;  -0.4;  +0.3;  +0.2];

%% ========================================================================
%  ETAPE 2 : CONSTRUCTION DE LA MATRICE D'ADMITTANCE Y_bus
%  ========================================================================
%  Regles (slides 18/70) :
%    - Diagonal Y_ii = somme de TOUTES les admittances connectees au noeud i
%                      (y compris les shunts vers la terre)
%    - Hors-diagonal Y_ij = -(admittance de la branche entre i et j)

yb = 1/Zb;  yc = 1/Zc;  yd = 1/Zd;
ye = 1/Ze;  yf = 1/Zf;
ya = 1/Za;  yg = 1/Zg;

Y = zeros(N, N);

% Hors-diagonaux
Y(3,2) = -yb;  Y(2,3) = Y(3,2);   % Branche b : bus 3 <-> bus 2
Y(3,1) = -yc;  Y(1,3) = Y(3,1);   % Branche c : bus 3 <-> bus 1
Y(2,1) = -yd;  Y(1,2) = Y(2,1);   % Branche d : bus 2 <-> bus 1
Y(2,4) = -ye;  Y(4,2) = Y(2,4);   % Branche e : bus 2 <-> bus 4
Y(1,4) = -yf;  Y(4,1) = Y(1,4);   % Branche f : bus 1 <-> bus 4

% Diagonaux (somme des admittances connectees, shunts inclus)
Y(1,1) = yc + yd + yf;         % Bus 1 : branches c, d, f
Y(2,2) = yb + yd + ye;         % Bus 2 : branches b, d, e
Y(3,3) = ya + yb + yc;         % Bus 3 : shunt a + branches b, c
Y(4,4) = ye + yf + yg;         % Bus 4 : branches e, f + shunt g

fprintf('=== Matrice d''admittance Y_bus (4x4) ===\n');
disp(Y);

%% ========================================================================
%  ETAPE 3 : S_complex - Formule compacte de puissance
%  ========================================================================
%  S = diag(V) * conj(Y * V)
%  car I = Y*V et S_j = V_j * conj(I_j)

S_complex = @(V_complex) diag(V_complex) * conj(Y * V_complex);

%% ========================================================================
%  ETAPE 4 : Conversion polaire -> complexe
%  ========================================================================
%  x = [v1; phi1; v2; phi2; v3; phi3; v4; phi4]
%  V_j = v_j * exp(j * phi_j)

polar_to_complex = @(x_vec) x_vec(1:2:end) .* exp(1i * x_vec(2:2:end));

%% ========================================================================
%  ETAPE 5 : Fonction de mismatch g(x) et bus PQ
%  ========================================================================

pq_buses = setdiff(1:N, slack_bus);  % = [2, 3, 4]
n_pq = length(pq_buses);             % = 3

g_func = @(x_vec) compute_mismatch(x_vec, Y, P_spec, Q_spec, ...
                                    pq_buses, S_complex, polar_to_complex);

%% ========================================================================
%  ETAPE 6 : Initialisation (condition "no load")
%  ========================================================================

x = zeros(2*N, 1);
for k = 1:N
    x(2*k-1) = 1.0;
    x(2*k)   = 0.0;
end

fprintf('=== Condition initiale (no load) ===\n');
for k = 1:N
    fprintf('  Bus %d : |V| = %.2f pu, angle = %.2f deg\n', ...
            k, x(2*k-1), rad2deg(x(2*k)));
end
fprintf('\n');

%% ========================================================================
%  ETAPE 7 : Boucle de Newton-Raphson
%  ========================================================================

tolerance = 1e-6;
max_iter = 50;
converged = false;
g_norm_history = [];

fprintf('=== Iterations de Newton-Raphson ===\n');
fprintf('%-6s %-15s\n', 'Iter', '||g(x)||_inf');
fprintf('-------------------------------\n');

for iter = 1:max_iter
    g_val = g_func(x);
    g_norm = norm(g_val, Inf);
    g_norm_history(end+1) = g_norm;
    fprintf('%-6d %-15.2e\n', iter, g_norm);
    
    if g_norm < tolerance
        converged = true;
        fprintf('\n*** Convergence en %d iterations ! ***\n\n', iter);
        break;
    end
    
    J = compute_jacobian_numerical(x, g_func, pq_buses, N);
    delta_x_red = J \ (-g_val);
    
    for k = 1:n_pq
        bus = pq_buses(k);
        x(2*bus-1) = x(2*bus-1) + delta_x_red(k);
        x(2*bus)   = x(2*bus)   + delta_x_red(n_pq + k);
    end
end

if ~converged
    fprintf('\n*** PAS DE CONVERGENCE apres %d iterations ***\n\n', max_iter);
end

%% ========================================================================
%  ETAPE 8 : Resultats
%  ========================================================================

v_mag   = x(1:2:end);
v_phase = x(2:2:end);
V_complex = polar_to_complex(x);

S_final = S_complex(V_complex);
P_calc = real(S_final);
Q_calc = imag(S_final);

fprintf('=== RESULTATS DU LOAD FLOW ===\n\n');
fprintf('%-6s %-12s %-14s %-12s %-12s %-8s\n', ...
        'Bus', '|V| (pu)', 'Angle (deg)', 'P (pu)', 'Q (pu)', 'Type');
fprintf('--------------------------------------------------------------\n');
for k = 1:N
    if k == slack_bus, btype = 'Slack'; else, btype = 'PQ'; end
    fprintf('%-6d %-12.6f %-14.6f %-12.6f %-12.6f  %s\n', ...
            k, v_mag(k), rad2deg(v_phase(k)), P_calc(k), Q_calc(k), btype);
end

fprintf('\nPuissance slack  : P = %.6f pu, Q = %.6f pu\n', ...
        P_calc(slack_bus), Q_calc(slack_bus));
fprintf('Pertes actives   : %.6f pu\n', sum(P_calc));
fprintf('Pertes reactives : %.6f pu\n', sum(Q_calc));

%% ========================================================================
%  ETAPE 9 : Courants et puissances dans les branches
%  ========================================================================

fprintf('\n=== COURANTS ET PUISSANCES DANS LES BRANCHES ===\n\n');

branches = {3, 2, yb, 'Yb (bus 3->2)';
            3, 1, yc, 'Yc (bus 3->1)';
            2, 1, yd, 'Yd (bus 2->1)';
            2, 4, ye, 'Ye (bus 2->4)';
            1, 4, yf, 'Yf (bus 1->4)'};

fprintf('%-18s %-14s %-14s %-14s\n', ...
        'Branche', '|I| (pu)', 'P_flow (pu)', 'Q_flow (pu)');
fprintf('------------------------------------------------------------\n');

for b = 1:size(branches, 1)
    i = branches{b, 1};
    j = branches{b, 2};
    y_br = branches{b, 3};
    nom = branches{b, 4};
    
    I_branch = y_br * (V_complex(i) - V_complex(j));
    S_branch = V_complex(i) * conj(I_branch);
    
    fprintf('%-18s %-14.6f %-14.6f %-14.6f\n', ...
            nom, abs(I_branch), real(S_branch), imag(S_branch));
end

%% ========================================================================
%  ETAPE 10 : Graphique de convergence
%  ========================================================================

figure;
semilogy(1:length(g_norm_history), g_norm_history, 'b-o', 'LineWidth', 2, ...
         'MarkerFaceColor', 'b');
grid on;
xlabel('Iteration');
ylabel('||g(x)||_{\infty}');
title('Convergence de Newton-Raphson - Reseau 4 bus');
hold on;
xl = xlim;
plot(xl, [tolerance tolerance], 'r--', 'LineWidth', 1.5);
text(xl(2), tolerance, ' Tolerance', 'VerticalAlignment', 'bottom', 'Color', 'r');
hold off;

%% ========================================================================
%  FONCTIONS AUXILIAIRES
%  ========================================================================

function g = compute_mismatch(x_vec, Y, P_spec, Q_spec, pq_buses, ...
                              S_complex_func, polar_to_complex_func)
    V_complex = polar_to_complex_func(x_vec);
    S = S_complex_func(V_complex);
    g_P = real(S(pq_buses)) - P_spec(pq_buses);
    g_Q = imag(S(pq_buses)) - Q_spec(pq_buses);
    g = [g_P; g_Q];
end

function J = compute_jacobian_numerical(x, g_func, pq_buses, N)
    n_pq = length(pq_buses);
    n_eq = 2 * n_pq;
    epsilon = 1e-7;
    idx_v   = 2*pq_buses - 1;
    idx_phi = 2*pq_buses;
    idx_all = [idx_v, idx_phi];
    n_var = length(idx_all);
    J = zeros(n_eq, n_var);
    for col = 1:n_var
        x_plus = x;
        x_plus(idx_all(col)) = x_plus(idx_all(col)) + epsilon;
        g_plus = g_func(x_plus);
        x_minus = x;
        x_minus(idx_all(col)) = x_minus(idx_all(col)) - epsilon;
        g_minus = g_func(x_minus);
        J(:, col) = (g_plus - g_minus) / (2 * epsilon);
    end
end