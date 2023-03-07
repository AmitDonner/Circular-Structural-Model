function [x_bar,y_bar,Sxx,Syy,Sxy]=GetMoments(x,y)
% Computing the moments of X and Y.
x_bar=mean(x);
y_bar=mean(y);
Sigma = cov(x,y);
Sxx = Sigma(1,1);
Syy = Sigma(2,2);
Sxy = Sigma(1,2);
end


