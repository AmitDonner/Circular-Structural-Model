function [lb1,ub1,lb2,ub2]=get_mu_bounds(Sxx,Syy,Sxy)
if Sxy<0 && Sxx-Syy<0 
    lb1=0*pi;
    ub1=0.25*pi;
    lb2=1*pi;
    ub2=1.25*pi;
    
elseif Sxy<0 && Sxx-Syy>0
    lb1=0.25*pi;
    ub1=0.5*pi;
    lb2=1.25*pi;
    ub2=1.5*pi; 
    
elseif Sxy>0 && Sxx-Syy>0
    lb1=0.5*pi;
    ub1=0.75*pi;
    lb2=1.5*pi;
    ub2=1.75*pi; 

else
    lb1=0.75*pi;
    ub1=1*pi;
    lb2=1.75*pi;
    ub2=2*pi; 
end
end