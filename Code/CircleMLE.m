function [sol,P_val,stat]=CircleMLE(x,y,r)
% checking input
if (all(size(x) ~= 1)) || (all(size(y) ~= 1))
    error('x and y must be numerical vectors')
end
% reshaping x and y to be column vectors.
x = reshape(x, [max(size(x)), 1]);
y = reshape(y, [max(size(y)), 1]);
% scale factor, if not supplied r is set to 1.
        if nargin < 3
            r = 1;
        end 
%% MLE using numerical solution
    % Steps:
    % 1. Computing moments of X,Y
    % 2. Testing for uniformity
        % 2.1. if Uniform - MLE for uniform case.
        % Else - Go to step 3. 
    % 3. Estimate the 2 possible values of mu.
    % 4. Defining R_1(kappa) and R_2(kappa) functions.
    % 5. Defining estimating equations of sigma^2 and rho 
    % 6. Estimating kappa for each possible value of mu. 
    % 7. Estimating  sigma, rho, a and b for each possible value of mu.
    % 8. Defining lower and upper bounds, and initial points for numerical solution
    % 9. Numerical solution
%% Step 1 - Computing moments of X,Y
[x_bar,y_bar,Sxx,Syy,Sxy]=GetMoments(x,y);

%% Step 2 - Testing for Uniformity
[P_val,stat] = UniformityTest(x,y,x_bar,y_bar,Sxx,Syy,Sxy);
if P_val > 0.05
    % Step 2.1
    sol = [0;0;CircleMLE_Uniform(x,y)];
else
%% Step 3 - Estimating mu using the Spectral decomposition 
    [Mu, lam]=FindMu(Sxx,Syy,Sxy);
    lambda_d = max(lam)-min(lam); 
%% Step 4 - R_1, R_2 functions.
    R1=@(kap) besseli(1,kap)./besseli(0,kap);
    R2=@(kap) besseli(2,kap)./besseli(0,kap);
    alpha=@(kap) 1-R1(kap).^2;
    beta=@(kap) R2(kap)-R1(kap).^2;
    
%% Step 5 - Estimating equations for rho and sigma^2.
    Rho=@(kap) sqrt( (max(lam)-min(lam)) ./ (R1(kap).^2-R2(kap)) );
    sigma_sqr=@(kap) 0.5*(Sxx + Syy) - 0.5*(1-R1(kap).^2) .* Rho(kap).^2; 

%% Step 6 - Estimating kappa for each value of mu.
    % moment estimator of centered MGF.
    
    M1=mean(exp(cos(Mu(1))*((x-x_bar)*r) + sin(Mu(1))*((y-y_bar)*r)));
    M2=mean(exp(cos(Mu(2))*((x-x_bar)*r) + sin(Mu(2))*((y-y_bar)*r)));
    % Estimating equation for kappa
    f1=@(k) exp(sigma_sqr(k)*r^2./2 - Rho(k)*r.*R1(k)) .* besseli(0,Rho(k)*r+k)./besseli(0,k)-M1;
    f2=@(k) exp(sigma_sqr(k)*r^2./2 - Rho(k)*r.*R1(k)) .* besseli(0,Rho(k)*r+k)./besseli(0,k)-M2;
    % Estimating kappa using numerical method
    opts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
                                   'MaxFunctionEvaluations', 1e10, ...
                                   'OptimalityTolerance', 1e-40, ...
                                   'StepTolerance', 1e-40, ...
                                   'Display','off');
    % Kappa lower bound 
    Kappa_lower_bound = @(kap) - alpha(kap)./beta(kap) - (Sxx + Syy)/lambda_d;
    kappa0 = bisection1(0, 500, 1e-10, Kappa_lower_bound, 100);
    ub = 500; lb = kappa0;
    try
        kappa_hat = [lsqnonlin(f1,kappa0,lb,ub,opts), lsqnonlin(f2,kappa0,lb,ub,opts)];
    catch 
        error('MME of kappa is incomputable, try using scaling factor r=1e-1, r=1e-2 or larger if needed.');
    end 
 %% 7. estimating  sigma, rho, a and b for each possible value of mu
    rho_hat = Rho(kappa_hat);
    sigma_hat=sigma_sqr(kappa_hat);
    if any(sigma_hat<0) 
        % if any sigma^2<0 then we will set it to be 1.
        sigma_hat(sigma_hat<0)=0;
    end
    % Convert to standard deviation 
    sigma_hat=sigma_hat.^(0.5);

    a_hat=x_bar-rho_hat.*cos(Mu).*R1(kappa_hat);
    b_hat=y_bar-rho_hat.*sin(Mu).*R1(kappa_hat);

%% Step 8 - Lower and upper bounds for optimization
    % lower and upper bounds for mu
    [l1, l2, u1, u2] = get_mu_bounds1(Sxy);
    % The order of the vector is: [kappa, mu, rho, a, b, sigma]

    lb1=[ 0, l1,   0, -inf, -inf,                 0];
    ub1=[ub, u1, inf,  inf,  inf, sqrt((Sxx+Syy)/2)];

    lb2=[ 0, l2,   0, -inf, -inf,                 0];
    ub2=[ub, u2, inf,  inf,  inf, sqrt((Sxx+Syy)/2)];

    % starting points (for each possible value of mu)
    x01=[kappa_hat(1); Mu(1); rho_hat(1); a_hat(1); b_hat(1); sigma_hat(1) + exprnd(1,1)];
    x02=[kappa_hat(2); Mu(2); rho_hat(2); a_hat(2); b_hat(2); sigma_hat(1) + exprnd(1,1)];

%% Step 9 - Minimizing the log likelihood numerically.
% -log likelihood function, VM is a function handle which we want to
    % minimize.
    n=size(x,1);
    VM = @(xx) VM_negLogLike3(xx,x,y,n); 
    % Optimization options
    options = optimoptions('fmincon','Display','off',...
                           'Algorithm','interior-point', ...
                           'OptimalityTolerance',1e-20, ...
                           'ConstraintTolerance',1e-20, ...
                           'MaxFunctionEvaluations',600);

    % minimization
    % Theta1,Theta2 are the solutions and fval1,fval2 are the -loglikelihood
    % values corresponding to Theta
    [Theta1,fval1]=fmincon(VM,x01,[],[],[],[],lb1,ub1,[],options);
    [Theta2,fval2]=fmincon(VM,x02,[],[],[],[],lb2,ub2,[],options);

    % Choosing optimal solution. 
    % The final solution (the vector sol) is the vector theta which minimize the 
    % -log likelihood value. 
    if fval2>fval1
        sol=reshape(Theta1, [6,1]);
    else
        sol=reshape(Theta2, [6,1]);
    end
end
end 
