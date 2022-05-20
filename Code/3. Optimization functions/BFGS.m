% BFGS.m
% "vanilla" BFGS optimization function
% function to find minimizer x of functional f

function [x,k,M] = BFGS(x,f,grad_f,M,tol) 
%{
INPUT:
    x       = initial guess for solution 
    f       = function handle of Objective Function f(x)
    grad_f  = function handle of Gradient of f(x)
    M       = inverse hessian approximation
    tol     = tolerance

OUTPUT: 
    x       = solution approximation
    k       = # algorithm iterations (for debugging)
    M       = inverse hessian approximation
%}
% Constants: 
    be = 1;             % contraction factor
    c1 = 0.125;         % Wolfe constant 
    
% Initial guess:
    x  = x(:);          % col vector
    fk = f(x);      
    gk = grad_f(x);
    pk = M*gk;          % search direction
    Tk = dot(gk,pk);    % tolerance parameter tau
    
% Begin:
    k  = 1;
    while sqrt(abs(Tk)) >= tol
        al  = 1.0;                  % init. step size: alpha = 1
        xup = x - al*pk;    
        
        % Backtracking Line Search: 
        while f(xup) > fk - c1*al*Tk               
            al  = al/(1 + be);      % step size alpha
            xup = x - al*pk;   
        end
        
        % Compute & Update: 
        sk = xup - x;
        yk = -gk;
        x  = xup;       
        gk = grad_f(x); 
        yk = yk + gk;
        rk = 1 / dot(yk,sk); 
        if k == 1
            M = (dot(yk,sk)/dot(yk,yk)) * eye(length(x)); % lambda*I
        end   
        M  = M - rk*(M*yk)*(sk');                         % + "right"
        M  = M - rk*sk*(((M')*yk)') + rk*sk*(sk');        % + "left" + 3rd
        fk = f(x);  
        pk = M*gk;
        Tk = dot(gk,pk); 
        k  = k + 1;
    end
    
end