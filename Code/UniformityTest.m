function [pval,stat]=UniformityTest(x,y,x_bar,y_bar,Sxx,Syy,Sxy)
n=length(x);
Var=mean(((x-x_bar).^2 + (y-y_bar).^2).^2)/2;
Stats=(Sxx-Syy)^2 + 4*Sxy^2;
stat=n*Stats/Var;
pval=chi2cdf(stat,2,'upper');
end