% LBFGS2.m
% L-BFGS optimization function we use with MPM
% function to find minimizer x of functional f

function [x,it,ng] = LBFGS2(x,f,grad_f,H,L2_rho,m,tol)
%{
INPUT:
    x       = initial guess for solution 
    f       = function handle of Objective Function f(x)
    grad_f  = function handle of Gradient of f(x)
    H       = function handle of inverse Hessian approximation
    L2_rho  = function handle to compute tolerance tau_MPM (MPM specific)
    m       = memory limit
    tol     = tolerance

OUTPUT: 
    x       = solution approximation
    it      = k = # algorithm iterations (for debugging)
    ng      = norm(grad_f) (for debugging)
%}

% Constants: 
    al = .25;           % Wolfe constant c1
    be = .5;            % contraction factor beta
    de = .9;
    
% Initialize:
    n  = size(x,1);     % assuming x is a col vector
    M  = 2*m;
    V  = zeros(n,M);    % "empty array" V
    S  = zeros(M,M);    % "empty array" Sigma
    
% Initial guess:
    fc      = f(x);
    gc      = grad_f(x);
    tau_MPM = L2_rho(x);
    p       = H(gc);        % H*grad_f = Id*grad_f = grad_f 
    tau     = dot(gc,p);    
    
% Begin:
    it      = 0;
    while( sqrt(abs(tau)) > tol && tau_MPM > n^-3.5)   % -3.5 best results
        dt      = 1;            % initial step size: alpha = 1
        xup     = x - dt*p;
        fup     = f(xup);
        gup     = grad_f(xup);
        tau_MPM = L2_rho(xup);
        
        % Backtracking Line Search: 
        while( fup > fc - al*dt*tau && abs(dot(gup,p)) > de*tau )
            dt  = dt/(1 + be);  % step size alpha
            xup = x - dt*p;
            fup = f(xup);
            gup = grad_f(xup);
        end
       
        % Compute & Update: 
        it  = it+1;
        s   = xup - x;
        y   = -gc;
        x   = xup;
        fc  = f(x);
        gc  = grad_f(x);
        y   = y + gc;
        tht = max( norm(y), norm(s) );   % avoiding dividing by 0
        y   = y / tht;
        s   = s / tht;
        r   = 1 / dot(y,s);
        
        % Redefine Function Handle: H = lambda*Id
        la  = dot(y,s) / dot(y,y);              % lambda  (scalar)
        H   = @(x) la*x;
        
        % Compute: 
        z   = H(y);
        V   = V - r*s*((y')*V);                 % *** edit ***
        W   = [z s];                            % new data (n x 2)
        L   = [0 -r; -r r*(1 + r*dot(z,y))];    % Lambda   (2 x 2)
                
        % Update: V & Sigma
        st  = min(2*(it-1) + 1,M+1);        % start index
        en  = min(2*it, M+2);               % end index
        V(:,st:en)      = W;                % concatenate new data
        S(st:en,st:en)  = L;
        V   = V(:, end-M+1:end);            % remove oldest data
        S   = S(end-M+1:end, end-M+1:end);
        
        % Update for next iteration:
        p   = (gc')*V;
        p   = S*(p');
        p   = V*p;
        p   = H(gc) + p;
        tau = dot(gc,p);
    end
    
    ng = max( abs(gc) );
end
