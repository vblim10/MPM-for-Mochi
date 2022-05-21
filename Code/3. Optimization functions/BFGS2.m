% BFGS2.m
% BFGS optimization function we use with MPM
% function to find minimizer x of functional f

function [x,k,M] = BFGS2(x,f,grad_f,M,L2_rho,tol) 
%{
INPUT:
    x       = initial guess for solution 
    f       = function handle of Objective Function f(x)
    grad_f  = function handle of Gradient of f(x)
    M       = inverse Hessian approximation (d*NG x d*NG)
    L2_rho  = function handle to compute tolerance tau_MPM (MPM specific)
    tol     = tolerance

OUTPUT: 
    x       = solution approximation
    k       = # algorithm iterations (for debugging)
    M       = inverse hessian approximation
%}
% Constants: 
    be = 0.50;              % contraction factor beta
    c1 = 0.125;             % Wolfe constant        
    
% Initial guess:
    x  = x(:);              % col vector
    fk = f(x);
    gk = grad_f(x);
    pk = M*gk;              % search direction
    Tk = dot(gk,pk);        % tolerance parameter tau
    tau_MPM = L2_rho(gk);           
    
% Begin: 
    k = 1;
    while tau_MPM > tol && sqrt(abs(Tk)) > tol
        al  = 1.0;              % init. step size: alpha = 1
        xup = x - al*pk;    
        
        % Backtracking Line Search:
        while f(xup) > fk - c1*al*Tk              
            al  = al/(1 + be);  % step size alpha
            xup = x - al*pk;   
        end
        
        % Compute & Update:
        sk = xup - x;
        yk = -gk;
        x  = xup;       
        gk = grad_f(x); 
        yk = yk + gk;
        rk = 1/dot(yk,sk); 
        if k == 1
            M = (dot(yk,sk)/dot(yk,yk))*eye(length(x)); % lambda*Id
        end   
        M  = M - rk*(M*yk)*(sk');                       % + "right"
        M  = M - rk*sk*(((M')*yk)') + rk*sk*(sk');      % + "left" + 3rd
        fk = f(x);  
        pk = M*gk;
        Tk = dot(gk,pk); 
        tau_MPM = L2_rho(gk); 
        k = k + 1;
    end

end