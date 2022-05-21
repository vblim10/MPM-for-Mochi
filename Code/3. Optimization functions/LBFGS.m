% LBFGS.m
% "vanilla" L-BFGS optimization function
% function to find minimizer x of functional f

function [x,k] = LBFGS(x,f,grad_f,H,m,tol) 
%{
INPUT:
    x       = initial guess for solution 
    f       = function handle of Objective Function f(x)
    grad_f  = function handle of Gradient of f(x)
    H       = function handle of inverse Hessian approx. (Identity)
    m       = memory limit
    tol     = tolerance

OUTPUT: 
    x       = solution approximation
    k       = # algorithm iterations (for debugging)
%}
% Constants: 
    be = 0.50;          % contraction factor beta
    c1 = 0.25;          % Wolfe constant 
    
% Initial guess:
    x  = x(:);          % col vector (n x 1)
    fk = f(x);      
    gk = grad_f(x);
    pk = H(gk);         % search direction
    Tk = dot(gk,pk);    % tolerance parameter tau
    
% Initialize:
    n  = size(x,1);     
    M  = 2*m;
    V  = zeros(n,M);    % "empty array" V
    S  = zeros(M,M);    % "empty array" Sigma
    
% Begin:
    k  = 1;
    while sqrt(abs(Tk)) > tol
        al  = 1.0;                  % init. step size: alpha = 1
        xup = x - al*pk;    
        
        % Backtracking Line Search: 
        while f(xup) > fk - c1*al*Tk               
            al  = al/(1 + be);      % step size alpha
            xup = x - al*pk;   
        end
        
        % Compute & Update: 
        k  = k + 1;
        sk = xup - x;
        yk = -gk;
        x  = xup;  
        fk = f(x);
        gk = grad_f(x); 
        yk = yk + gk;
        rk = 1 / dot(yk,sk); 
        
        % Redefine Function Handle: H = lambda*Id
        la = dot(yk,sk) / dot(yk,yk);               % lambda  (scalar)
        H  = @(x) la*x;
        
        % Compute: 
        zk = H(yk);
        W  = [zk sk];                               % new data (n x 2)
        La = [0 -rk; -rk rk*(1 + rk*dot(zk,yk))];   % Lambda   (2 x 2)
        
        % Update: V & Sigma
        st = min(2*(k-1) + 1,M+1);          % start index
        en = min(2*k, M+2);                 % end index
        V(:,st:en)     = W;                 % concatenate new data
        S(st:en,st:en) = La;
        V  = V(:, end-M+1:end);             % remove oldest data
        S  = S(end-M+1:end, end-M+1:end);
        
        % Update for next iteration:
        pk = (gk')*V;
        pk = S*(pk');
        pk = V*pk;
        pk = H(gk) + pk;
        Tk = dot(gk,pk);          
    end
    
end