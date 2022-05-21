% nH_skewColl_2d.m 
% 2d neo-Hookean model (skew impact of cyllinders with gravity)

d = 2;              % dimension
CFL = 0.25;
tf  = 12.0;         % final time
dtk = 0.10;         % uniform time step
g = 9.80665/100;    % gravitational constant

% Material & Plot settings
ppc  = 2;           % # particles per grid cell = ppc^d
psize = 25;         % for plotting particles
sqSz = 1;           % for plotting grid nodes

% Elastic constants
% (E = 31.685, nu = 0.44022) -> mu = 11.0001, Ka = 88.3378, La = 81.0044
E_ym = 31.685; nu_pr = 0.44022;
mu = 0.5*E_ym/(1+nu_pr);    % shear modulus
Ka = (E_ym/3)/(1-2*nu_pr);  % bulk modulus 
La = Ka-(2/3)*mu;           % Lame's 1st parameter
md_p = 5.0;                 % particle mass density (for all p)
    
% Elastic params struct for funct. handles
elast.mu = mu; elast.Ka = Ka; elast.La = La;    % struct for funct. handles

% ======================== EULERIAN MESH =================================
a = 20.0; b = 12.0; % Eul mesh domain: x=0:a, y=0:b
h = [1.0, 1.0];     % grid stepsize (20x12 grid cells)   
% h = [0.5, 0.5];     % grid stepsize (40x24 grid cells)

% Grid initialization
x = 0:h(1):a;       % without padding
y = 0:h(2):b;
nx = length(x)-1;   % # x-axis partitions
ny = length(y)-1;   % # y-axis partitions

[Ei,Ej] = meshgrid(x,y);    % grid node positions
NG = length(Ei(:));         % # grid nodes

% Eulerian grid struct
grid.O = [min(x) min(y)]; % origin = grid bottom-left corner (1x2)
grid.h = h;             % 1x2
grid.n = [nx, ny];      % 1x2
grid.height = Ej(:);    % NGx1

% Boundary Condition grid indices
BCg = zeros(size(Ei,1),size(Ei,2));
BCg(1,:)=1; BCg(end,:)=1; BCg(:,1)=1; BCg(:,end)=1;
BCg = find(BCg(:));

% ======================= LAGRANGIAN PARTICLES ===========================
% Cyllinder initial values
    % c1 
    c1x = 3.0;   c1y = 5.0;     % center = (3.0,5.0)
    v1x = 0.75;  v1y = 0.0;     % velocity = <0.75,0>
    r1in = 1.0;  r1out = 2.0;   % inner, outer radius
    % c2 (*** assuming: c2 & c1 have same radii ***)
    c2x = 16.0;  c2y = 7.0;     % center = (16.0,7.0)
    v2x = -0.75; v2y = 0.0;     % velocity = <-0.75,0>
    r2in = 1.0;  r2out = 2.0;   % inner, outer radius

% c1 (center=(3.0,5.0), radius=2.0)
    xc = (c1x-r1out) + h(1)/(2*ppc) : h(1)/ppc : (c1x+r1out) - h(1)/(2*ppc);  
    yc = (c1y-r1out) + h(2)/(2*ppc) : h(2)/ppc : (c1y+r1out) - h(2)/(2*ppc);
    [Xi1,Xj1] = meshgrid(xc,yc);  % c1: particle positions (rect. box)

    % Remove particles "outside" of c1 (center=(3.0,5.0), radius=2.0)
    outside = ( Xi1(:)-c1x ).^2 + ( Xj1(:)-c1y ).^2;
    outside = [find( double(outside > r1out^2) );% outside particle indices
               find( double(outside < r1in^2) )];
    Xi1 = Xi1(:); Xj1 = Xj1(:);                  % reshape to col vectors
    Xi1(outside) = []; Xj1(outside) = [];        % remove outside particles
    
% c2 is c1 shifted (*** assuming: c2 & c1 have same radius ***)
    Xi2 = Xi1 + (c2x-c1x);
    Xj2 = Xj1 + (c2y-c1y);
    
% All particles coords 
    Phi_k = [Xi1, Xj1; Xi2, Xj2];   % NPx2
    NP = size(Phi_k,1);             % # particles
    
% =========================== CHECK GRAPH ================================
    figure(1)
    scatter(Phi_k(:,1),Phi_k(:,2),psize,'filled','c'); hold on;
    scatter(Ei(:),Ej(:),'+','b'); hold off;

% =========================== Step k = 0 =================================
tic
k = 0;
t = 0;      % time at step k

% Compute Particle Volumes (Tot_Volume/NP)
    Vol = prod(h)*(1/ppc^d)*ones(NP,1); 
    
% Set Particle Mass (md_p = 5.0 for all p)
    mass_p = md_p*ones(NP,1).*Vol;  
    
% Initial Values (k=0)
    % Deformation gradient (F = Id for all particles) 
    Id = eye(d);
    Fk = repmat(Id(:),[1,NP]);  % 4xNP
    Jk = det_2d(Fk);            % 1xNP
    
    % APIC var
    Bm = 0*Fk;                  % 4xNP
    
    % Particle Velocity (2xNP)
    %         c1                 c2
    Uk = [v1x*ones(1,NP/2), v2x*ones(1,NP/2); 
          v1y*ones(1,NP/2), v2y*ones(1,NP/2)];     
    
    % "Actual" Energy over time 
    % funct. handle: Vol, mass_p, elast, & g are constant over time
    RealEngy = @(J,Fk,Uk,Xp) RealEngy_comp(J,Fk,Uk,Xp,Vol,mass_p,elast,g); 
    
    % inverse Hessian approx. (for BFGS)
%     Hinv = eye(d*NG);
       
    % Tolerance (for minimization)
    tol = min([max(h)^4, 1e-4, max(h)^-4]); % cases: h<1, h=1, h>1
    
    % ANIMATION
    T = [t,t]; 
    PHX = [Phi_k(:,1), Phi_k(:,1)]; 
    PHY = [Phi_k(:,2), Phi_k(:,2)]; 
    JPK = [det_2d(Fk)', det_2d(Fk)'];
    
% ============================= Step k ===================================
while t < tf
    % "Actual" Energy over time 
    EvT(k+1,:) = [t,RealEngy(Jk,Fk,Uk,Phi_k)]; % [t,E,Q,K,El,G] 
    
    %  APIC CFL (non-uniform time steps)
%     dtk = CFL*h(1)/( max(sum(abs(Uk),1) + (6/h(1))*sqrt(d)*sqrt(sum(Bm.^2,1))) );
%     dtk = min(dtk,tf-t);
    
    % APIC TRANSFER particle to grid 
    mom_p = bsxfun(@times,Uk,mass_p');      % particle mom    2xNP
    [mom_g,mass_g] = P2G_APIC(grid,Phi_k,Bm,mom_p,mass_p);
    mass_g = mass_g + eps;  % so we're not dividing mom_g by zeros 
    u = bsxfun(@rdivide, mom_g, mass_g');   % vel = mom/mass  2xNG
    u(:,BCg) = 0;                           % enforce BC
    
    % DYNAMICS 
    % Initial guess (velocity)     
    w = u(:);           % (d*NG)x1
    
    % Redefine Function Handles  
    Hinv0 = @(w) w;     % init inv Hessian approx = Id (for LBFGS)
    E = @(w) E_comp(w,u,dtk,Fk,Vol,mass_g,grid,Phi_k,elast,g); 
    gradE = @(w) gradE_comp(w,u,dtk,Fk,Vol,mass_g,grid,Phi_k,elast,g);
    L2_rho = @(gradE) L2_rho_comp(gradE,mass_g,h);
    
    % BFGS: [x,k,M] = BFGS2(x,f,grad_f,M,L2_rho,tol) 
%     [w,kb,Hinv] = BFGS2(w,E,gradE,Hinv,L2_rho,tol);

    % LBFGS: [x,it,ng] = LBFGS2(x,f,grad_f,H,L2_rho,m,tol)
    [w,kb,ng] = LBFGS2(w,E,gradE,Hinv0,L2_rho,5,tol); % m=5 (memory limit)

    % APIC TRANSFER grid to particle (advection) 
    w = reshape(w,[d NG]);                      % 2xNG
    w(:,BCg) = 0;                               % enforce BC
    xup = [Ei(:)'; Ej(:)'] + dtk*( (w+u)/2 );   % 2xNG
    [Phi_up,Vk,Bm] = G2P_APIC(grid,Phi_k,xup,w); 

    % GRAPH 
    scatter(Phi_up(:,1),Phi_up(:,2),psize,'filled','co'); hold on;
    scatter(Ei(:),Ej(:),'+','b'); hold on;
    title(num2str(kb),num2str(t)); xlabel(num2str(k)); ylabel(num2str(dtk));
    hold off; pause(1e-10);
    
    % ANIMATION
    T = [T,t]; 
    PHX = [PHX, Phi_up(:,1)]; 
    PHY = [PHY, Phi_up(:,2)]; 
    JPK = [JPK, det_2d(Fk)'];
    
    % Update vars for next time step 
    k = k + 1;
    t = t + dtk;
    Fk = F_comp(grid,Phi_k,dtk,u+w,Fk); % 4xNP
    Jk = det_2d(Fk);                    % 1xNP
    Phi_k = Phi_up;                     % NPx2
    Uk = Vk;                            % 2xNP
    
end
time = toc;

% "Actual" Energy over time 
EvT(k+1,:) = [t,RealEngy(Jk,Fk,Uk,Phi_k)]; % [t,E,Q,K,El,G] 

% ++++++++++++++++++++++++++++ PLOTS +++++++++++++++++++++++++++++++++++++

% Energy Plots
figure(2)
plot(EvT(:,1),EvT(:,2),'r');   hold on; % Total Energy
plot(EvT(:,1),EvT(:,3),'-go'); hold on; % Potential Energy
plot(EvT(:,1),EvT(:,4),'-bs'); hold on; % Kinetic Energy
legend('Total Energy','Potential Energy','Kinetic Energy');
title('Energy');
xlabel('time t'); hold off;

figure(3)
plot(EvT(:,1),EvT(:,2),'r');    % Total Energy
title('Total Energy'); 

% ------------------------- ANIMATION DATA -------------------------------
% Write Energy data to text file ([t,E,Q,K,El,G])
writematrix(EvT,'Engy.dat','Delimiter',';');  % (k+1 x 6)

% Write Particle positions to text file
writematrix(PHX,'PHX.dat','Delimiter',';');  % (x-coords) (NP x k+1)
writematrix(PHY,'PHY.dat','Delimiter',';');  % (y-coords) (NP x k+1)

% Write J=det(Fp) particle data to text file
writematrix(JPK,'Jpk.dat','Delimiter',';');  % (NP x k+1)

% Write Grid Node positions to text file
writematrix([Ei(:),Ej(:)],'GNxy.dat','Delimiter',';');  % (NG x 2)

% -------------------------- ANIMATED SIM --------------------------------
% SAVE Workspace and IMPORT into nH_skewColl_2d_ani.m  
% in order to render simulation as gif/mp4

% ++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++

function [Engy] = RealEngy_comp(J,Fk,Uk,Xp,Vol,mass_p,elast,g)
% function to compute "actual" energy at time tk
%{  
    input var dims:
    J           1 x NP
    Fk          4 x NP
    Uk          2 x NP
    Xp          NP x 2
    Vol         NP x 1
    mass_p      NP x 1
    elast       struct of scalar constants

    funct var dims: 
    C               4 x NP
    trC             1 x NP
    H               1 x NP
    E,Q,K,El,G      scalar
    Engy            1 x 5
%}
% dimension
    d = 2;

% Elastic constants
    mu = elast.mu; % shear modulus
    La = elast.La; % Lame's 1st parameter

% Compute:          
    C = Tmult_2d(Fk,Fk); 
    trC = C(1,:) + C(4,:);
    
% Potential Energy 
    % H(F) = W(F'F) = W(C)
    H = 0.5*mu*(trC - d) - mu*log(J) + 0.5*La*log(J).*log(J);% JST17 & mpm course eqn. 46
    El = dot(H,Vol');           % Elastic potential energy 
    G = g*dot(mass_p,Xp(:,2));  % Gravitational potential energy
    Q = El + G;                 % Potential Energy

% Kinetic Energy
    K = 0.5*(dot(Uk,Uk,1) * (mass_p));
    
% Total Energy 
    E = Q + K; 
    
% Output
    Engy = [E, Q, K, El, G];
end

function [E,El,K,G] = E_comp(w,u,dtk,Fk,Vol,mass_g,grid,Xp,elast,g)
% function to compute E, given w
%{  
    input var dims:
    w,u         2 x NG
    Fk          4 x NP
    Vol         NP x 1
    mass_g      NG x 1

    funct var dims: 
    Ftk         4 x NP
    J           1 x NP
    C           4 x NP
    trC         1 x NP
    H           1 x NP
    E,K,El,G    scalar
%}
% Elastic constants
    mu = elast.mu; % shear modulus
    La = elast.La; % Lame's 1st parameter

    d = 2;
    NG = length(mass_g);
    w = reshape(w,[d NG]);
    u = reshape(u,[d NG]);
    
% Compute Ftk
    Ftk = F_comp(grid,Xp,dtk/2,w+u,Fk);
    
% Elastic Energy: El
    % H(F) = W(F'F) = W(C)
    J = det_2d(Ftk);
    C = Tmult_2d(Ftk,Ftk);
    trC = C(1,:) + C(4,:);
    H = 0.5*mu*(trC - 2) - mu*log(J) + 0.5*La*log(J).*log(J);% JST17 & mpm course eqn. 46
    El = dot(H,Vol');    
    
% Gravitational Energy: G
    G = g*dot(mass_g,grid.height);
    
% Kinetic Energy: K
    K = (dot(w-u,w-u,1) * mass_g)/8;
    
% Total Energy: E = El + 0.25*dt*G + K
    E = El + 0.25*dt*G + K;
end

function [gradE] = gradE_comp(w,u,dtk,Fk,Vol,mass_g,grid,Xp,elast,g)
% function to compute gradE, given w
%{  
    input var dims:
    w,u         2 x NG
    Fk          4 x NP
    Vol         NP x 1
    mass_g      NG x 1

    funct var dims: 
    Ftk,P,B,FinvT               4 x NP
    J                           1 x NP
    gradEl,gradG,gradK,gradE    (2*NG) x 1
%}
% Elastic constants
    mu = elast.mu; % shear modulus
    La = elast.La; % Lame's 1st parameter

    d = 2;
    NG = length(mass_g);
    w = reshape(w,[d NG]);
    u = reshape(u,[d NG]);
    
% Compute Ftk
    Ftk = F_comp(grid,Xp,dtk/2,w+u,Fk);
    
% Elastic Energy: gradEl
    % P(F) = DH(F)
    J = det_2d(Ftk);
    FinvT = invT_2d(Ftk);
    % JST17 & mpm course, eqn. 48
    P = mu*(Ftk - FinvT) + bsxfun(@times, La*log(J), FinvT);
    
    % B
    B = multT_2d(P,Fk); % 4xNP  
    B = 0.25 * dtk * bsxfun(@times,B,Vol');
    
    % [vg] = P2G_vec(grid,Xp,vp,b)
    gradEl = P2G_vec(grid,Xp,B(1:2,:),[1 0]) + P2G_vec(grid,Xp,B(3:4,:),[0 1]); % dxNG
    gradEl = gradEl(:);   % (d*NG)x1
    
% Gravitational Energy: gradG
    gradG = [zeros(1,NG); g*(mass_g')]; % positive b/c it is on the RHS of vel. update
    gradG = gradG(:);
    
% Kinetic Energy: gradK
    gradK = 0.25*bsxfun(@times,(w-u),(mass_g'));  % dxNG
    gradK = gradK(:);     % (d*NG)x1
    
% Total Energy: gradE
    gradE = gradEl + 0.25*dt*gradG + gradK;
end

function [tau] = L2_rho_comp(gradE,mass_g,h)
    % mass_g     NGx1
    d = 2;
    NG = length(mass_g);
    gradE = reshape(gradE,[d NG]); % dxNG
    tau = sum(gradE.^2,1).*(mass_g')/prod(h);
    tau = sqrt(sum(tau));
end

function [detA] = det_2d(A)
% A     4 x n
% detA  1 x n
detA = A(1,:).*A(4,:) - A(2,:).*A(3,:);
end

function [AinvT] = invT_2d(A)
% A         4 x n
% AinvT     4 x n
% detA      1 x n
detA = A(1,:).*A(4,:) - A(2,:).*A(3,:);
detA = detA + eps;  % avoid dividing by zero

AinvT = [ A(4,:); 
         -A(3,:);
         -A(2,:);
          A(1,:)];
AinvT = bsxfun(@rdivide,AinvT,detA);
end

% vectorized (and FASTER than C):
function [M] = mult_2d(A,B)
% M,A,B:     4 x n

% M = A * B
M = [A(1,:).*B(1,:) + A(3,:).*B(2,:);
     A(2,:).*B(1,:) + A(4,:).*B(2,:);
     A(1,:).*B(3,:) + A(3,:).*B(4,:);
     A(2,:).*B(3,:) + A(4,:).*B(4,:)];
end

function [M] = multT_2d(A,B)
% M,A,B:     4 x n

% M = A * B'
M = [A(1,:).*B(1,:) + A(3,:).*B(3,:);
     A(2,:).*B(1,:) + A(4,:).*B(3,:);
     A(1,:).*B(2,:) + A(3,:).*B(4,:);
     A(2,:).*B(2,:) + A(4,:).*B(4,:)];
end

function [M] = Tmult_2d(A,B)
% M,A,B:     4 x n

% M = A' * B
M = [A(1,:).*B(1,:) + A(2,:).*B(2,:);
     A(3,:).*B(1,:) + A(4,:).*B(2,:);
     A(1,:).*B(3,:) + A(2,:).*B(4,:);
     A(3,:).*B(3,:) + A(4,:).*B(4,:)];
end

% ========================= Wrapper Functions =============================
function fp = G2P_PIC(grid,Xp,fg,b)
%{  
var dimensions:
    Xp:     NP x 2
    fg:     1 x NG
    b =     [bx,by]
    fp:     1 x NP  
%}

    % Parse inputs (unpack)
    Ox = grid.O(1); Oy = grid.O(2);
    hx = grid.h(1); hy = grid.h(2);
    nx = grid.n(1); ny = grid.n(2);
    Yp = Xp(:,2); Xp = Xp(:,1);
    bx = b(1); by = b(2);

    % Call C function
    fp = G2P_PIC_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, fg, bx,by);
end

function fg = P2G_PIC(grid,Xp,fp,b)
%{  
var dimensions:
    Xp:     NP x 2
    fp:     1 x NP
    b =     [bx,by]

    fg:     1 x NG  
%}

    % Parse inputs (unpack)
    Ox = grid.O(1); Oy = grid.O(2);
    hx = grid.h(1); hy = grid.h(2);
    nx = grid.n(1); ny = grid.n(2);
    Yp = Xp(:,2);   Xp = Xp(:,1);
    bx = b(1);      by = b(2);

    % Call C function
    fg = P2G_PIC_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, fp, bx,by);
end

function [vp] = G2P_vec(grid,Xp,vg,b)
%{  
var dimensions:
    Xp:     NP x 2
    vg:     2 x NG
    b =     [bx,by]

    vp:     2 x NP  
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2);
hx = grid.h(1); hy = grid.h(2);
nx = grid.n(1); ny = grid.n(2);
Yp = Xp(:,2);   Xp = Xp(:,1);
fg1 = vg(1,:);  fg2 = vg(2,:);
bx = b(1);      by = b(2);

% Call C function
[fp1,fp2] = G2P_vec_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, fg1,fg2, bx,by);

% Package Outputs
vp = [fp1; fp2];
end

function [vg] = P2G_vec(grid,Xp,vp,b)
%{  
var dimensions:
    Xp:     NP x 2
    vp:     2 x NP
    b =     [bx,by]

    vg:     2 x NG  
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2);
hx = grid.h(1); hy = grid.h(2);
nx = grid.n(1); ny = grid.n(2);
Yp = Xp(:,2);   Xp = Xp(:,1);
fp1 = vp(1,:);  fp2 = vp(2,:);
bx = b(1);      by = b(2);

% Call C function
[fg1,fg2] = P2G_vec_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, fp1,fp2, bx,by);

% Package Outputs
vg = [fg1; fg2];
end

function [Fup] = F_comp(grid,Xp,dt,v,F)
%{  
var dimensions:
    Xp:     NP x 2
    dt:     scalar
    v:      2 x NG
    F:      4 x NP

    Fup:    4 x NP
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2);
hx = grid.h(1); hy = grid.h(2);
nx = grid.n(1); ny = grid.n(2);
Yp = Xp(:,2);   Xp = Xp(:,1);
vx = v(1,:);    vy = v(2,:);
F1 = F(1,:);    F2 = F(2,:);
F3 = F(3,:);    F4 = F(4,:);

% Call C function
[Fup1,Fup2,Fup3,Fup4] = F_comp_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, dt, vx,vy, F1,F2,F3,F4);

% Package Outputs
Fup = [Fup1; Fup2; Fup3; Fup4];
end

function [Phi,Vp,B] = G2P_APIC(grid,Xp,xup,w)
%{  
var dimensions:
    Xp:     NP x 2
    xup,w:  2 x NG
    
    Phi:    NP x 2
    Vp:     2 x NP
    B:      4 x NP
%}

% Parse inputs (unpack)
Ox = grid.O(1);     Oy = grid.O(2);
hx = grid.h(1);     hy = grid.h(2);
nx = grid.n(1);     ny = grid.n(2);
Yp = Xp(:,2);       Xp = Xp(:,1);
xup1 = xup(1,:);    xup2 = xup(2,:);
w1 = w(1,:);        w2 = w(2,:);

% Call C function
[Phi1,Phi2, Vp1,Vp2, B1,B2,B3,B4] = G2P_APIC_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, xup1,xup2, w1,w2);

% Package Outputs
Phi = [Phi1, Phi2];
Vp = [Vp1; Vp2];
B = [B1; B2; B3; B4];
end

function [mom_g,mass_g] = P2G_APIC(grid,Xp,B,mom_p,mass_p)
%{  
var dimensions:
    Xp:         NP x 2
    B:          4 x NP
    mom_p:      2 x NP
    mass_p:     NP x 1
    
    mom_g:      2 x NG
    mass_g:     NG x 1
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2);
hx = grid.h(1); hy = grid.h(2);
nx = grid.n(1); ny = grid.n(2);
Yp = Xp(:,2);   Xp = Xp(:,1);
B1 = B(1,:);    B2 = B(2,:);
B3 = B(3,:);    B4 = B(4,:);
mom_p1 = mom_p(1,:); mom_p2 = mom_p(2,:);
mass_p = mass_p(:);

% Call C function
[mom_g1,mom_g2, mass_g] = P2G_APIC_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, B1,B2,B3,B4, mom_p1,mom_p2, mass_p);

% Package Outputs
mom_g = [mom_g1; mom_g2];
end
