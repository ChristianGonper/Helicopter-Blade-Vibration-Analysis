function [M, K] = assembleMK(gamma, m_sym, EI_sym, G_sym, R)
% assembleMK: Assembles symbolic M and K matrices
% gamma: vector of symbolical functions of x
% m_sym, EI_sym, G_sym: symbolic functions/values

syms x;
N = length(gamma);  % Number of assumed functions
M = sym(zeros(N, N));  % Mass matrix
K = sym(zeros(N, N));  % Stiffness matrix

for i = 1:N
    for j = 1:N
        M(i,j) = int(m_sym * gamma(i) * gamma(j), x, 0, 1) * R;
        term1 = (EI_sym / R^3) * diff(gamma(i), x, 2) * diff(gamma(j), x, 2);
        term2 = (G_sym / R) * diff(gamma(i), x, 1) * diff(gamma(j), x, 1);
        K(i,j) = int(term1 + term2, x, 0, 1);
    end
end