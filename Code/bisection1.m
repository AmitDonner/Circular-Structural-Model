function [res,flag]=bisection1(a,b,tol,f,max_iter)
    % f is a function handle
    % a - left side, b - right side
    % tol = tolerance
    % max_iter - maximum possible iteration 
    %% Initial points check
    if f(a)*f(b)>0 
        while sign(f(a))==sign(f(b)) && a <= b
            a=a+0.01;
        end       
    end
    if f(a)*f(b)>0
        res=1;
        flag=1;
    else
        %% Iterations
        flag=0;
        iter=1; 
        % 3 arrays documenting the values of a,b,c
        c=mean([a,b]);
        while iter<=max_iter && abs(f(c))>=tol
            iter=iter+1;
            if sign(f(c))==sign(f(a)) 
                a=c;
            else
                b=c;
            end
            c=mean([a,b]);
        end
        if iter>max_iter
            flag = 1; 
        end 
        res=c;
    end 
  
end




