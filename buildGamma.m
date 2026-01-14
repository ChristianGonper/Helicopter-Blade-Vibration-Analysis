function gamma = buildGamma(N, bcType)
% buildGamma: Generates N assumed functions for Rayleigh-Ritz
% satisfying all BCs specified in the lecture notes

syms x a b;
gamma = sym(zeros(1, N));

if strcmpi(bcType, 'hinged')

    % Mode 1: Rigid body mode satisfies all 4 CCs
    gamma(1) = x;

    % Higher modes: x^n + a*x^(n+1) + b*x^(n+2)
    % Starting from n=3 ensures S(0)=0 and S''(0)=0
    for i = 2:N
        n = i + 1;
        f = x^n + a*x^(n+1) + b*x^(n+2);

        % Natural BCs at tip (x=1)
        d2f = diff(f, x, 2);
        d3f = diff(f, x, 3);

        eqs = [subs(d2f, x, 1) == 0, subs(d3f, x, 1) == 0];
        sol = solve(eqs, [a, b]);

        gamma(i) = simplify(subs(f, [a, b], [sol.a, sol.b]));
    end

elseif strcmpi(bcType, 'cantilevered')

    % Starting from n=2 ensures S(0)=0 and S'(0)=0
    for i = 1:N
        n = i + 1;
        f = x^n + a*x^(n+1) + b*x^(n+2);

        % Natural BCs at tip (x=1)
        d2f = diff(f, x, 2);
        d3f = diff(f, x, 3);

        eqs = [subs(d2f, x, 1) == 0, subs(d3f, x, 1) == 0];
        sol = solve(eqs, [a, b]);

        gamma(i) = simplify(subs(f, [a, b], [sol.a, sol.b]));
    end
else
    error('Unknown BC type');
end