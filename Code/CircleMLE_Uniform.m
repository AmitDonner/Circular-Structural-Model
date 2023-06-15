function sol = CircleMLE_Uniform(x,y)
a0 = mean(x);
b0 = mean(y);
rho0 = mean(sqrt((x-a0).^2 + (y-b0).^2));
sigma0 = sqrt(0.5 * mean((x-a0).^2 + (y-b0).^2 - rho0^2));
x0 = [rho0; a0; b0; sqrt(sigma0)];
%% Step 2 - Minimizing the log likelihood numerically.
% -log likelihood function, VM is a function handle which we want to
    % minimize.
n=size(x,1);
VM=@(xx)VM_negLogLikeUniform(xx,x,y,n); 
% Optimization options
options = optimoptions('fmincon','Display','off',...
    'Algorithm','interior-point');

% minimization
% Lower and upper bounds
lb1 = [0  , -inf, -inf,   0];
ub1 = [inf,  inf,  inf, inf];
sol = fmincon(VM,x0,[],[],[],[],lb1,ub1,[],options);
end 