function [y, f] = VM_negLogLike3(par,x,y,n)
%  kappa= par(1); mu= par(2); rho=par(3); a=par(4); b=par(5); sigma=par(6) (no sigma^2)
kappa = par(1); 
mu = par(2); 
rho = par(3); 
a = par(4); 
b = par(5); 
sigma = par(6);
% x - the x sample
% y - the y sample
%% D = sqrt(A^2+B^2)
A = kappa.* cos(mu) + (x - a).*rho./sigma.^2;
B = kappa.* sin(mu) + (y - b).*rho./sigma.^2;
D = sqrt(A.^2 + B.^2);
%% negative log-likelihood with exponential scaling and correction
% besseli(nu,z,1) for exponential scaling
y = -sum(log(besseli(0,D,1))) - sum(D) + n*log(besseli(0,kappa,1)) + n*kappa + ...
      (1/(2*sigma^2)) * sum((x-a).^2 + (y-b).^2 + rho^2) + ...
      n * log(sigma^2);
%% Gradient
if nargout > 1 % gradient required
    r = 1/2 * (1 - besseli(2, D)./besseli(0, D));
    
    f = zeros([6,1]); % column vector to be filled with the gradient's components.
    % dl/dkappa
    f(1) = -sum(r .* (kappa + rho./sigma.^2 .* (cos(mu).*(x-a) + sin(mu).*(y-b))) ...
                - besseli(1, kappa)./besseli(0, kappa));
    % dl/dmu
    f(2) = -sum(r .* kappa.*rho./sigma.^2 .* (-sin(mu).*(x-a) + cos(mu).*(y-b)));
    % dl/drho
    f(3) = -sum(r .* (kappa./sigma.^2 .*(cos(mu).*(x-a) + sin(mu).*(y-b)) + ...
                (rho./sigma.^4) .* ((x-a).^2 + (y-b).^2)) ...
                -rho/sigma.^2);
    % dl/da
    f(4) = -sum(-r.* ((kappa.*rho.*cos(mu))./sigma.^2+ (rho.^2 * (x - a))./sigma.^4)  + (x-a)./sigma.^2);
    % dl/db
    f(5) = -sum(-r.* ((kappa.*rho.*sin(mu))./sigma.^2+ (rho.^2 * (y - b))./sigma.^4)  + (y-b)./sigma.^2);
    % dl/dsigma^2
    f(6) = -sum(-r .* (kappa*rho./sigma.^4 .* (cos(mu).*(x - a) + sin(mu).*(y - b)) ...
            + rho.^2/sigma.^6 .* ((x-a).^2 + (y-b).^2)) ...
            + 0.5*sigma^(-4) .* ((x-a).^2 + (y-b).^2 +rho^2) - 1./sigma.^2);
end


