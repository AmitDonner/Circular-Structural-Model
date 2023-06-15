function [mu,lambda]=FindMu(Sxx, Syy, Sxy)
    %% Spectral Decompose
    Sigma = [Sxx, Sxy; Sxy, Syy];
    [V,D]=eig(Sigma,'vector');
    [~,ind]=min(D);
    %% Finding mu
    mu_temp(1)=mod(atan(V(2,ind)/V(1,ind)) + 2*pi, 2*pi);
    mu_temp(2)=mod(mu_temp(1)+pi,pi*2);
    mu = sort(mu_temp);
    lambda = D;
end 