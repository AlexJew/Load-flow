function J = make_jacobian_matrix(V_abs, V_angle, Y_bus, n_nodes)
% Compute the jacobian matrix for a grid to facilitate the computation of a
% load flow problem.
% 
% - IMPORTANT - this function assumes that the slack bus is node 1!
%
%
% INPUTS
%  V_abs (Nx1 matrix, with N number of nodes, including slack bus)
%   Nodal voltage magnitude (pu)
%  V_angle (Nx1 matrix)
%   Phase angle (rad) of nodal voltages
%  Y_bus (NxN matrix)
%   Bus admittance matrix in per unit.
%  n_nodes (scalar value)
%   Number of nodes in the grid (including the slack bus)!

% OUTPUT
% J ( 2(N-1) x 2(N-1) matrix )
%   Jacobian matrix. This matrix contains the ``sensitivity coefficients'' for
%   all buses of the networks except for the slack bus. The structure of
%   the matrix is as follows:
%   
%   --                                                                         --
%   | dP2/dv2     dP2/dv3     ...     dP2/dvN     dP2/dPhi1   ...     dPN/dPhiN |
%   | dP3/dv2     dP3/dv3     ...     dP3/dvN     dP3/dPhi1   ...     dPN/dPhiN |
%          .....
%   | dPN/dv2     dPN/dv3     ...     dPN/dvN     dPN/dPhi1   ...     dPN/dPhiN |
%   | dQ2/dv2     dQ2/dv3     ...     dQ2/dvN     dQ2/dPhi1   ...     dQ2/dPhiN |
%   | dQ3/dv2     dQ3/dv3     ...     dQ3/dvN     dQ3/dPhi1   ...     dQ3/dPhiN | 
%          .....
%   | dQN/dv2     dQN/dv3     ...     dQN/dvN     dQN/dPhi1   ...     dQN/dPhiN | 
%   --                                                                         --
% where Pj and Qj are active and reactive nodal power injections at node j,
% vj is the magnitude of the voltage phasor at node j, and PhiN is the
% angle of the voltage phasor at node j.
%
% Fabrizio Sossan (HESSO VS), 2023


Yabs = abs(Y_bus);
Yangle = angle(Y_bus);

dPdV = zeros(n_nodes,n_nodes);
dPdPhi = zeros(n_nodes,n_nodes);
dQdV = zeros(n_nodes,n_nodes);
dQdPhi = zeros(n_nodes,n_nodes);

for i = 1:n_nodes

    for h = 1:n_nodes

        if i ~= h

            dPdV(i,h) = Yabs(i,h)*V_abs(i)*cos(V_angle(i)-V_angle(h)-Yangle(i,h));
            % Error below?
            %dPdPhi(i,h) = Yabs(i,h)*V_abs(i)*V_abs(h)*sin(V_angle(i)-V_angle(h)-Yangle(i,h));
            % Correction:
            dPdPhi(i,h) = -Yabs(i,h)*V_abs(i)*V_abs(h)*sin(V_angle(i)-V_angle(h)+Yangle(i,h));

			dQdV(i,h) = Yabs(i,h)*V_abs(i)*sin(V_angle(i)-V_angle(h)+Yangle(i,h));
            dQdPhi(i,h) = -Yabs(i,h)*V_abs(i)*V_abs(h)*cos(V_angle(i)-V_angle(h)-Yangle(i,h));

        elseif i == h
            sum = 0; sumT=0; sumQ = 0; sumQT = 0;

            for j = 1:n_nodes

                sum = sum + Yabs(i,j)*V_abs(j)*cos(V_angle(i)-V_angle(j)-Yangle(i,j));
                sumT = sumT +(-Yabs(i,j)*V_abs(i)*V_abs(j)*sin(V_angle(i)-V_angle(j)-Yangle(i,j)));

                sumQ = sumQ + Yabs(i,j)*V_abs(j)*sin(V_angle(i)-V_angle(j)-Yangle(i,j));
                sumQT = sumQT +(Yabs(i,j)*V_abs(i)*V_abs(j)*cos(V_angle(i)-V_angle(j)-Yangle(i,j)));

            end

            dPdV(i,i) = Yabs(i,i) *V_abs(i)*cos(Yangle(i,i)) +  sum;
            dPdPhi(i,i) = sumT - (-V_abs(i)^2*Yabs(i,i)*sin(-Yangle(i,i)));

            dQdV(i,i) = -Yabs(i,i) *V_abs(i)*sin(Yangle(i,i)) +  sumQ;
            dQdPhi(i,i) = sumQT - (V_abs(i)^2*Yabs(i,i)*cos(-Yangle(i,i)));
        end
    end

end

J = [dPdV(2:n_nodes, 2:n_nodes) dPdPhi(2:n_nodes, 2:n_nodes); dQdV(2:n_nodes, 2:n_nodes) dQdPhi(2:n_nodes, 2:n_nodes)];

end