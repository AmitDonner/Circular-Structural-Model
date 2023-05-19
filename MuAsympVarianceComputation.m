syms t1 t2 sigma rho kappa mu % Defining symbolic variables
% Assumptions regarding symbolic variables
assume(t1, 'real'); % t1, t1 are real numbers
assume(t2, 'real');
assume(sigma>0); % sigma,kappa,rho are positive 
assume(rho>0);
assume(kappa > 0);
assume(0<=mu & mu<=pi); % mu \in [0,pi].
% R_k(kappa) symbolic functions, k=1,2 
R1 = besseli(1, kappa)/besseli(0, kappa);
R2 = besseli(2, kappa)/besseli(0, kappa);

% Symbolic centered moments generating function
MomFun = exp(sigma^2/2 * (t1^2 + t2^2) - rho * R1 *(t1 *cos(mu) + t2 * sin(mu))) * ...
    (besseli(0, sqrt(rho^2*(t1^2 + t2^2) +2*rho*kappa *(t1 *cos(mu) + t2 * sin(mu))+ kappa^2))/besseli(0, kappa)); 

% Symbolic diffrentiation of MomFun
cMOM = @(dx, dy) symfun(diff(diff(MomFun, t1, dx), t2, dy), [t1 t2]);
 
% cMOM(k,j) = E(x-mu_x)^k * (y- mu_y)^j
E_x4 = cMOM(4,0); % E(x-mu_x)^4
E_y4 = cMOM(0,4); % E(y-mqu_y)^4
sig_x = cMOM(2,0); % E(x-mu_x)^2
sig_y = cMOM(0,2); % E(y-mu_y)^2
sig_xy = cMOM(1,1); % Cov(x,y)
E_x2y2 = cMOM(2,2); % E(x-mu_x)^2(y-mu_y)^2
E_x3y = cMOM(3,1); % E(x-mu_x)^3(y-mu_y)
E_xy3 = cMOM(1,3); % E(x-mu_x)(y-mu_y)^3

% Building the asymptotic covariance
% matrix, \boldsymbol{C} of Sxx, Syy, Sxy.
D11 = simplify(expand(E_x4(0,0) - (sig_x(0,0))^2), 'IgnoreAnalyticConstraints', true);
D21 = simplify(expand(E_x2y2(0,0)-sig_x(0,0)*sig_y(0,0)), 'IgnoreAnalyticConstraints', true);
D22 = simplify(expand(E_y4(0,0)-(sig_y(0,0))^2), 'IgnoreAnalyticConstraints', true);
D31 = simplify(expand(E_x3y(0,0)-sig_xy(0,0)*sig_x(0,0)), 'IgnoreAnalyticConstraints', true);
D32 = simplify(expand(E_xy3(0,0) - sig_xy(0,0)*sig_y(0,0)), 'IgnoreAnalyticConstraints', true);
D33 = simplify(expand(E_x2y2(0,0) - (sig_xy(0,0))^2), 'IgnoreAnalyticConstraints', true);

% The matrix \boldsymbol{C}
Cov = [D11, D21, D31;
       D21, D22, D32;
       D31, D32, D33];

% Vector of partial derivative of \Phi(\mu) w.r.t Var(X), Var(Y), Cov(X,Y).
Nabla_g = [-sin(2*mu);
           sin(2*mu);
           2*cos(2*mu)];

% \sigma_{\mu}^2 comutation
VarTheta = (Nabla_g') * Cov * Nabla_g; % asymptotic variance of the numerator in eq (5.15)
VarTheta = simplify(simplify(expand(simplify(expand(VarTheta),...
           'steps', 10)),'Steps',50), 'Steps',40); % simplification 
VarDenom = (-2 * rho^2 * (R2 - R1^2))^2;
% \sigma_\mu^2 computation
VarMu = simplify(VarTheta/VarDenom, 'IgnoreAnalyticConstraints', true, 'Steps', 80);
