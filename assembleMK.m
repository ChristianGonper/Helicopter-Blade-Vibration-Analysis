function [M, K] = assembleMK(gamma, m_sym, EI_sym, G_sym, R)
% assembleMK: Assembles symbolic M and K matrices
% gamma: vector of symbolical functions of x
% m_sym, EI_sym, G_sym: symbolic functions/values

syms x;
N = length(gamma);  % Number of assumed functions
M = sym(zeros(N, N));  % Mass matrix
K = sym(zeros(N, N));  % Stiffness matrix

% Following energy method formulas:
% Potential Energy V = 1/2 * int(EI * (S'')^2 + G * (S')^2, r, 0, R)
% Kinetic Energy T = 1/2 * int(m * (w_dot)^2, r, 0, R)

% Adimensionalizing r = x*R, dr = R*dx
% V = 1/2 * int(EI(x) * (d2S/dx2 * 1/R^2)^2 + G(x) * (dS/dx * 1/R)^2, x, 0, 1) * R
% V = 1/2 * int(EI(x)/(R^3) * (S'')^2 + G(x)/R * (S')^2, x, 0, 1)
% T = 1/2 * R * int(m(x) * omega^2 * S^2, x, 0, 1)

for i = 1:N
    for j = 1:N
        M(i,j) = int(m_sym * gamma(i) * gamma(j), x, 0, 1) * R;
        term1 = (EI_sym / R^3) * diff(gamma(i), x, 2) * diff(gamma(j), x, 2);
        term2 = (G_sym / R) * diff(gamma(i), x, 1) * diff(gamma(j), x, 1);
        K(i,j) = int(term1 + term2, x, 0, 1);
    end
end