% CHANNEL_FLOW  Newtonian vs shear-thinning (power-law) plane channel flow.
% MATLAB port of src/channelflow/solver.py (finite volume + Picard iteration).
% NOTE: written to mirror the validated Python code; run it once in your own
% MATLAB install and compare against results/table1_fixed_flowrate.csv.
%
%   d/dy[ mu(|u'|) u' ] = -G,   u'(0)=0 (symmetry),  u(h)=0
clear; clc; close all;

h = 1; K = 1; Ub = 1;  N = 200;  nlist = [1 0.8 0.6 0.4 0.2];
figure; hold on; grid on;
for n = nlist
    G  = K*(Ub*(2*n+1)/(n*h^((n+1)/n)))^n;          % G for a given bulk velocity
    mu = @(gd) K*(gd.^2 + 1e-16).^((n-1)/2);         % regularised power law
    [y,u,it] = solve_channel(N,G,mu,h,1e-10,20000,0.5);
    p  = (n+1)/n;
    ue = n/(n+1)*(G/K)^(1/n)*(h^p - y.^p);           % analytical
    err = norm(u-ue)/norm(ue);
    fprintf('n=%.1f  G=%.4f  iters=%d  L2 error=%.2e\n',n,G,it,err);
    plot(u/trapz(y,u)*h, y, 'LineWidth',1.8,'DisplayName',sprintf('n = %.1f',n));
end
xlabel('u/U_b'); ylabel('y/h'); legend show; title('Velocity profiles at equal flow rate');

function [y,u,it] = solve_channel(N,G,mu,h,tol,maxit,relax)
    y  = linspace(0,h,N+1)'; dy = y(2)-y(1); M = N;
    u  = G/(2*mu(1))*(h^2 - y.^2);
    for it = 1:maxit
        mf = mu(abs(diff(u))/dy);                    % face viscosities (N x 1)
        main = zeros(M,1); up = zeros(M-1,1); lo = zeros(M-1,1);
        main(1) = -2*mf(1);  up(1) = 2*mf(1);
        main(2:M) = -(mf(2:M) + mf(1:M-1));
        up(2:M-1) = mf(2:M-1);
        lo(1:M-1) = mf(1:M-1);
        A   = spdiags([[lo;0], main, [0;up]], -1:1, M, M);
        rhs = -G*dy^2*ones(M,1);
        un  = [A\rhs; 0];
        unext = (1-relax)*u + relax*un;
        res = max(abs(unext-u))/max(abs(unext));
        u = unext;
        if res < tol, break; end
    end
end
