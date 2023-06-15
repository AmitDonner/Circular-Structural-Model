function [lb1, lb2, ub1, ub2] = get_mu_bounds1(Sxy)
if Sxy<0  
    lb1=0;
    ub1=0.5*pi;
    lb2=pi;
    ub2=1.5*pi;
else
    lb1=0.5*pi;
    ub1=pi;
    lb2=1.5*pi;
    ub2=2*pi; 

end