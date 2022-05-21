% VE_skewColl_3d.m 
% 3d Viscoelastic model (skew impact of balls with gravity)

d = 3;              % dimension
CFL = 0.25;
tf  = 25.0;         % final time
dtk = 0.10;         % uniform time step
g = 9.80665/100;    % gravitational constant

% Material & Plot settings
ppc  = 2;           % # particles per grid cell = ppc^d
psize = 25;         % for plotting particles
sqSz = 1;           % for plotting grid nodes

% Viscoelastic constants
    % Hyperelastic (rubber)
    %{
    E_ym = 31.685; nu_pr = 0.44022;
    mu = 0.5*E_ym/(1+nu_pr);    % shear modulus
    Ka = (E_ym/3)/(1-2*nu_pr);  % bulk modulus 
    La = Ka-(2/3)*mu;           % Lame's 1st parameter
    mu_N = 0;                   % viscous Newton constant
    Wi = 5000;                  % Weissenberg number
    md_p = 5.0;                 % particle mass density (for all p)
    %}
    
    % Shaving foam
    %
    mu = 5;             % shear modulus
    La = 50;            % Lame's 1st parameter
    mu_N = 1.0*10^-4;   % viscous Newton constant
    Wi = 0.5;           % Weissenberg number
    md_p = 0.2;         % particle mass density (for all p)
    %}
    
    % Toothpaste
    %{
    mu = 0.839;         % shear modulus
    La = 8.39;          % Lame's 1st parameter
    mu_N = 1.0*10^-1;   % viscous Newton constant
    Wi = 0.4;           % Weissenberg number
    md_p = 1.0;         % particle mass density (for all p)
    %}
    
    % Viennetta ice cream
    %{
    mu = 1;             % shear modulus
    La = 10;            % Lame's 1st parameter
    mu_N = 5.0*10^-5;   % viscous Newton constant
    Wi = 0.1;           % Weissenberg number
    md_p = 1.0;         % particle mass density (for all p)
    %}  

% VE params struct for funct. handles
Qconst.mu = mu;     Qconst.La = La; 
Qconst.mu_N = mu_N; Qconst.Wi = Wi;

% ======================== EULERIAN MESH =================================
a = 20.0; b = 12.0; c = 8.0; % Eul mesh domain: x=0:a, y=0:b, z=0:c
h = [1.0, 1.0, 1.0];         % grid stepsize (20x12x8 grid cells)   
% h = [0.5, 0.5, 0.5];         % grid stepsize (40x24x16 grid cells) 

% Grid initialization
x = 0:h(1):a;       % without padding
y = 0:h(2):b;
z = 0:h(3):c;
nx = length(x)-1;   % # x-axis partitions
ny = length(y)-1;   % # y-axis partitions
nz = length(z)-1;   % # z-axis partitions

[Ei,Ej,Ek] = meshgrid(x,y,z);   % grid node positions
NG = length(Ei(:));             % # grid nodes
    
% Eulerian grid struct
grid.O = [min(x) min(y) min(z)]; % origin = grid bottom-left corner (1x3)
grid.h = h;                      % 1x3
grid.n = [nx, ny, nz];           % 1x3
grid.height = Ek(:);             % NGx1

% Boundary Condition grid indices
BCg = zeros(size(Ei,1),size(Ei,2),size(Ei,3));
BCg(1,:,:)=1; BCg(end,:,:)=1; 
BCg(:,1,:)=1; BCg(:,end,:)=1;
BCg(:,:,1)=1; BCg(:,:,end)=1;
BCg = find(BCg(:));
    
% ======================= LAGRANGIAN PARTICLES ===========================
% Ball initial values
    % c1 
    c1x = 3.0;   c1y = 5.0;   c1z = 4.0;    % center = (3.0,5.0,4.0)
    v1x = 0.75;  v1y = 0.0;   v1z = 0.5;    % velocity = <0.75,0,0.5>
    r1in = 0.0;  r1out = 2.0;               % inner, outer radius
    % c2 (*** assuming: c2 & c1 have same radii ***)
    c2x = 15.0;  c2y = 7.0;   c2z = 4.0;    % center = (15.0,7.0,4.0)
    v2x = -0.75; v2y = 0.0;   v2z = 0.5;    % velocity = <-0.75,0,0.5>
    r2in = 0.0;  r2out = 2.0;               % inner, outer radius
    
% c1 (center=(3.0,5.0,4.0), radius=2.0)
    xc = (c1x-r1out) + h(1)/(2*ppc) : h(1)/ppc : (c1x+r1out) - h(1)/(2*ppc);  
    yc = (c1y-r1out) + h(2)/(2*ppc) : h(2)/ppc : (c1y+r1out) - h(2)/(2*ppc);
    zc = (c1z-r1out) + h(3)/(2*ppc) : h(3)/ppc : (c1z+r1out) - h(3)/(2*ppc);
    [Xi1,Xj1,Xk1] = meshgrid(xc,yc,zc); % c1 particle positions (rect. box)
    
    % Remove particles "outside" of c1 (center=(3.0,5.0,4.0), radius=2.0)
    outside = ( Xi1(:)-c1x ).^2 + ( Xj1(:)-c1y ).^2 + ( Xk1(:)-c1z ).^2;
    outside = [find( double(outside > r1out^2) );% outside particle indices
               find( double(outside < r1in^2) )];
    Xi1 = Xi1(:); Xj1 = Xj1(:); Xk1 = Xk1(:);    % reshape to col vectors
    Xi1(outside)=[]; Xj1(outside)=[]; Xk1(outside)=[];   % remove particles
    
% c2 (just c1 shifted) (*** assuming: c2 & c1 have same radius ***)
    Xi2 = Xi1 + (c2x-c1x);
    Xj2 = Xj1 + (c2y-c1y);
    Xk2 = Xk1 + (c2z-c1z);
    
% All particles coords 
    Phi_k = [Xi1, Xj1, Xk1; Xi2, Xj2, Xk2]; % NPx3           
    NP = size(Phi_k,1);                     % # particles
    
% =========================== CHECK GRAPH ================================
    figure(1)
    scatter3(Phi_k(:,1),Phi_k(:,2),Phi_k(:,3),psize,'filled','c'); hold on;
    scatter3(Ei(:),Ej(:),Ek(:),sqSz,'+','b'); hold off;

% ============================ Step k = 0 ================================
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
    Fk = repmat(Id(:),[1,NP]);      % 9xNP
    Jk = det_3d(Fk);                % 1xNP
    
    % Elastic Stress vars
    Bob_k = repmat(Id(:),[1,NP]);   % 9xNP
    Bk = B_comp(Fk,Bob_k);          % 9xNP
    
    % APIC var
    Bm = 0*Fk;                      % 9xNP
    
    % Particle Velocity (3xNP)
    %         c1                 c2
    Uk = [v1x*ones(1,NP/2), v2x*ones(1,NP/2); 
          v1y*ones(1,NP/2), v2y*ones(1,NP/2);
          v1z*ones(1,NP/2), v2z*ones(1,NP/2)];     
    
    % "Actual" Energy over time 
    % funct. handle: grid, Vol, mass_p, & Qconst are constant over time
    RealEngy = @(J,B_El,Uk,w,Xp) RealEngy_comp(J,B_El,Uk,w,Xp,grid,Vol,mass_p,Qconst,g);
    u = P2G_vec(grid,Phi_k,Uk,[0 0 0]); w = u; % so we can compute init. Engy
       
    % Tolerance (for minimization)
    tol = min([max(h)^4, 1e-4, max(h)^-4]); % cases: h<1, h=1, h>1
    
    % ANIMATION
    T = [t,t]; 
    PHX = [Phi_k(:,1), Phi_k(:,1)]; 
    PHY = [Phi_k(:,2), Phi_k(:,2)]; 
    PHZ = [Phi_k(:,3), Phi_k(:,3)]; 
    JPK = [det_3d(Fk)', det_3d(Fk)'];

% ============================= Step k ===================================
while t < tf
    % "Actual" Energy over time 
    EvT(k+1,:) = [t,RealEngy(Jk,Bk,Uk,w,Phi_k)]; % [t,E,Q,K,Q_El,Q_N,G]
    
    %  APIC CFL (non-uniform time steps)
%     dtk = CFL*h(1)/( max(sum(abs(Uk),1) + (6/h(1))*sqrt(d)*sqrt(sum(Bm.^2,1))) );
%     dtk = min(dtk,tf-t);
    
    % APIC TRANSFER particle to grid 
    mom_p = bsxfun(@times,Uk,mass_p');      % particle mom    3xNP
    [mom_g,mass_g] = P2G_APIC(grid,Phi_k,Bm,mom_p,mass_p);
    mass_g = mass_g + eps;  % so we're not dividing mom_g by zeros 
    u = bsxfun(@rdivide, mom_g, mass_g');   % vel = mom/mass  3xNG
    u(:,BCg) = 0;                           % enforce BC
    
    % DYNAMICS
    % Initial guess (velocity)     
    w = u(:);           % (d*NG)x1
    
    % Redefine Function Handles  
    Hinv0 = @(w) w;     % init inv Hessian approx = Id (for LBFGS)
    E = @(w) E_comp(w,u,dtk,Fk,Bob_k,Vol,mass_g,grid,Phi_k,Qconst,g); 
    gradE = @(w) gradE_comp(w,u,dtk,Fk,Bob_k,Vol,mass_g,grid,Phi_k,Qconst,g);
    L2_rho = @(gradE) L2_rho_comp(gradE,mass_g,h);

    % LBFGS: [x,it,ng] = LBFGS2(x,f,grad_f,H,L2_rho,m,tol)
    [w,kb,ng] = LBFGS2(w,E,gradE,Hinv0,L2_rho,5,tol); % m=5 (memory limit)

    % APIC TRANSFER grid to particle (advection) 
    w = reshape(w,[d NG]);                  % 3xNG
    w(:,BCg) = 0;                           % enforce BC
    xup = [Ei(:)'; Ej(:)'; Ek(:)'] + dtk*w;
    [Phi_up,Vk,Bm] = G2P_APIC(grid,Phi_k,xup,w); 

    % GRAPH 
    scatter3(Phi_up(:,1),Phi_up(:,2),Phi_up(:,3),psize,'filled','co'); hold on; % particles
    scatter3(Ei(:),Ej(:),Ek(:),sqSz,'+','b'); hold on; % grid
    title(num2str(kb),num2str(t)); xlabel(num2str(k)); ylabel(num2str(dtk));
    hold off; pause(1e-10);
    
    % ANIMATION
    T = [T,t]; 
    PHX = [PHX, Phi_up(:,1)]; 
    PHY = [PHY, Phi_up(:,2)]; 
    PHZ = [PHZ, Phi_up(:,3)];
    JPK = [JPK, det_3d(Fk)'];
    
    % Update vars for next time step 
    k = k + 1;
    t = t + dtk;
    Fk = F_comp(grid,Phi_k,dtk,w,Fk);   % 9xNP
    Jk = det_3d(Fk);                    % 1xNP
    Bob_k = Bob_comp(Bob_k,dtk,Qconst,grid,Phi_k,w); % 9xNP
    Bk = B_comp(Fk,Bob_k);              % 9xNP
    Phi_k = Phi_up;                     % NPx3
    Uk = Vk;                            % 3xNP

end
time = toc;

% "Actual" Energy over time 
EvT(k+1,:) = [t,RealEngy(Jk,Bk,Uk,w,Phi_k)]; % [t,E,Q,K,Q_El,Q_N,G]

% ++++++++++++++++++++++++++++ PLOTS +++++++++++++++++++++++++++++++++++++

% Energy Plots
figure(2)
plot(EvT(:,1),EvT(:,2),'r');   hold on; % Total Energy
plot(EvT(:,1),EvT(:,3),'-go'); hold on; % Potential Energy
plot(EvT(:,1),EvT(:,4),'-bs'); hold on; % Kinetic Energy
legend('Total Energy','Potential Energy','Kinetic Energy');
xlabel('time t'); hold off;
title('Energy');

figure(3)
plot(EvT(:,1),EvT(:,2),'r');    % Total Energy
title('Total Energy'); 

% ------------------------- ANIMATION DATA -------------------------------
% Write Energy data to text file ([t,E,Q,K,Q_El,Q_N,G])
writematrix(EvT,'Engy.dat','Delimiter',';');  % (k+1 x 7)

% Write Particle positions to text file
writematrix(PHX,'PHX.dat','Delimiter',';');  % (x-coords) (NP x k+1)
writematrix(PHY,'PHY.dat','Delimiter',';');  % (y-coords) (NP x k+1)
writematrix(PHZ,'PHZ.dat','Delimiter',';');  % (z-coords) (NP x k+1)

% Write J=det(Fp) particle data to text file
writematrix(JPK,'Jpk.dat','Delimiter',';');  % (NP x k+1)

% Write Grid Node positions to text file
writematrix([Ei(:),Ej(:),Ek(:)],'GNxyz.dat','Delimiter',';');  % (NG x 3)

% -------------------------- ANIMATED SIM --------------------------------
% SAVE Workspace and IMPORT into VE_skewColl_3d_ani.m  
% in order to render simulation as gif/mp4

% ++++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++

function [Engy] = RealEngy_comp(J,Bk,Uk,w,Xp,grid,Vol,mass_p,Qconst,g)
% function to compute "actual" energy at time tk
%{  
    input var dims:
    J           1 x NP
    B_El        9 x NP
    Uk          3 x NP
    w           3 x NG
    Vol         NP x 1
    mass_p      NP x 1
    Qconst      struct of scalar constants

    funct var dims: 
    Vol_k               NP x 1
    trB_El              1 x NP
    H                   1 x NP
    DV,S                9 x NP
    E,Q,K,Q_El,Q_N      scalar
    Engy                1 x 5
%}
% dimension
    d = 3;

% Elastic & Viscous constants
    mu = Qconst.mu;     % shear modulus
    La = Qconst.La;     % Lame's 1st parameter
    mu_N = Qconst.mu_N; % viscous Newton constant
    
% Potential Energy 
    Vol_k = (J').*Vol;
    trB = Bk(1,:) + Bk(5,:) + Bk(9,:);
    % H(B_El)
    H = 0.5*mu*(trB - d) - mu*log(J) + 0.5*La*(J-1).*(J-1);
    Q_El = dot(H,Vol');     % elastoplastic portion of potential energy
    DV = DV_comp(grid,Xp,w);
    S = 0.5*(DV + tran_3d(DV));
    Q_N = mu_N*( dot(S,S,1) *Vol_k ); % Newtonian viscous potential energy
    G = g*dot(mass_p,Xp(:,3));        % Gravitational potential energy
    Q = Q_El + Q_N + G;               % Potential Energy
    
% Kinetic Energy
    K = 0.5*(dot(Uk,Uk,1) * (mass_p));
    
% Total Energy 
    E = Q + K; 
    
% Output
    Engy = [E, Q, K, Q_El, Q_N, G];
end

function [E,El,En,K] = E_comp(w,u,dt,Fk,Bob_k,Vol,mass_g,grid,Xp,Qconst,g)
% function to compute E, given w
%{  
    input var dims:
    w,u         3 x NG
    Fk,Bob_k    9 x NP
    Vol         NP x 1
    mass_g      NG x 1

    funct var dims: 
    Dw          9 x NP
    Id          9 x NP
    Fup         9 x NP
    Jup         1 x NP
    Bob_up      9 x NP
    detBob      1 x NP
    dk          1 x NP
    Bup         9 x NP
    trB         1 x NP
    PB          1 x NP
    Dw_sym      9 x NP
    R           1 x NP
    El,En,K,E   scalar
%}
% Elastic & Viscous constants
    mu = Qconst.mu;     % shear modulus
    La = Qconst.La;     % Lame's 1st parameter
    Wi = Qconst.Wi;     % Weissenberg number
    mu_N = Qconst.mu_N; % viscous Newton constant
    
    d = 3;
    NG = length(mass_g);
    NP = length(Vol);
    w = reshape(w,[d NG]);
    u = reshape(u,[d NG]);
    
% Compute
    Dw = DV_comp(grid,Xp,w);
    Id = repmat([1 0 0 0 1 0 0 0 1]',1,NP);
    Fup = mult_3d(Id + dt*Dw ,Fk);
    Jup = det_3d(Fup);
    Bob_up = Bob_k + dt*( mult_3d(Dw,Bob_k) + multT_3d(Bob_k,Dw) + (1/Wi)*(Id-Bob_k) );
    detBob = det_3d(Bob_up) + eps;          % 1xNP
    dk = ( (Jup.*Jup) ./ detBob ).^(1/3);   % 1xNP
    Bup = bsxfun(@times,dk,Bob_up);
    
% Elastic Energy: El
    trB = Bup(1,:) + Bup(5,:) + Bup(9,:);                % 1xNP
    PB = 0.5*mu*(trB - d) - mu*log(Jup) + La*(Jup-1).^2; % 1xNP
    El = dot(PB,Vol');   
    
% Viscous Energy: En
    Dw_sym = 0.5*( Dw + tran_3d(Dw) );  % 9xNP
    R = sum(Dw_sym.*Dw_sym,1);          % 1xNP
    Vol_up = (Jup').*Vol;               % NPx1
    En = mu_N * dot(R,Vol_up');

% Gravitational Energy: G
    G = g*dot(mass_g,grid.height);
    
% Kinetic Energy: K
    K = (dot(w-u,w-u,1) * mass_g)/2;
    
% Total Energy: E = El + dt*En + dt*G + K
    E = El + dt*En + dt*G + K;
end

function [gradE] = gradE_comp(w,u,dt,Fk,Bob_k,Vol,mass_g,grid,Xp,Qconst,g)
% function to compute gradE, given w
%{  
    input var dims:
    w,u         3 x NG
    Fk,Bob_k    9 x NP
    Vol         NP x 1
    mass_g      NG x 1

    funct var dims: 
    Dw          9 x NP
    Id          9 x NP
    Fup         9 x NP
    Jup         1 x NP
    Bob_up      9 x NP
    detBob      1 x NP
    dk          1 x NP
    Bup         9 x NP
    P,P1,P2     9 x NP
    trP         1 x NP
    Bob_up_invT 9 x NP
    Dw_sym      9 x NP
    Q           9 x NP
    gradPot, gradK, gradE   (3*NG) x 1
%}
% Elastic & Viscous constants
    mu = Qconst.mu;     % shear modulus
    La = Qconst.La;     % Lame's 1st parameter
    Wi = Qconst.Wi;     % Weissenberg number
    mu_N = Qconst.mu_N; % viscous Newton constant

    d = 3;
    NG = length(mass_g);
    NP = length(Vol);
    w = reshape(w,[d NG]);
    u = reshape(u,[d NG]);
    
% Compute
    Dw = DV_comp(grid,Xp,w);
    Id = repmat([1 0 0 0 1 0 0 0 1]',1,NP);
    Fup = mult_3d(Id + dt*Dw ,Fk);
    Jup = det_3d(Fup);
    Bob_up = Bob_k + dt*( mult_3d(Dw,Bob_k) + multT_3d(Bob_k,Dw) + (1/Wi)*(Id-Bob_k) );
    detBob = det_3d(Bob_up) + eps;          % 1xNP
    dk = ( (Jup.*Jup) ./ detBob ).^(1/3);   % 1xNP
    Bup = bsxfun(@times,dk,Bob_up);
    
% Elastic Energy: gradEl
    P = 0.5*mu*Bup+ 0.5*bsxfun(@times,( La*(Jup.*Jup - Jup)-mu ),Id); % 9xNP
    % gradE_OB_1
    trP = P(1,:) + P(5,:) + P(9,:);             % 1xNP
    Bob_up_invT = invT_3d(Bob_up);              % 9xNP
    P1 = 2*multT_3d(invT_3d(Fup),Fk) - multT_3d(Bob_up_invT,Bob_k) - Tmult_3d(Bob_k,Bob_up_invT);
    P1 = (dt/3) * bsxfun(@times,trP,P1);        % 9xNP
    % gradE_OB_2
    P2 = 2*dt*P; 
    
% Viscous Energy: gradEn
    Dw_sym = 0.5*( Dw + tran_3d(Dw) );          % 9xNP
    Q = 2*mu_N*Dw_sym;                          % 9xNP
    
% Potential Energy: gradPot = gradEl + gradEn 
    % interpolate: R = (P1 + P2 + dt*Q) * Vol'
    R = bsxfun(@times, P1 + P2 + dt*Q, Vol');   % 9xNP
    
    % [vg] = P2G_vec(grid,Xp,vp,b)
    gradPot = P2G_vec(grid,Xp,R(1:3,:),[1 0 0]) + P2G_vec(grid,Xp,R(4:6,:),[0 1 0]) + P2G_vec(grid,Xp,R(7:9,:),[0 0 1]); % 3xNG
    gradPot = gradPot(:);   % (3*NG)x1

% Gravitational Energy: gradG
    gradG = dt*[zeros(1,NG); zeros(1,NG); g*(mass_g')]; % positive b/c it is on the RHS of vel. update
    gradG = gradG(:);
    
% Kinetic Energy: gradK
    gradK = bsxfun(@times,(w-u),(mass_g'));  % dxNG
    gradK = gradK(:);     % (d*NG)x1
    
% Total Energy: gradE
    gradE = gradPot + gradG + gradK;
end

function [tau] = L2_rho_comp(gradE,mass_g,h)
    % mass_g     NGx1 
    d = 3;
    NG = length(mass_g);
    gradE = reshape(gradE,[d NG]);              % dxNG
    tau = sum(gradE.^2,1).*(mass_g')/prod(h);   % 1xNG
    tau = sqrt(sum(tau));                       % scalar
end

% Viscoelastic:
function [Bob_up] = Bob_comp(Bob_k,dt,Qconst,grid,Xp,w)
%{  
    input var dims:
    Bob_k       9 x NP
    dt          scalar
    Qconst      struct of scalar constants
	grid        struct for interpolation
    Xp          NP x 3
    u           3 x NG
    w           3 x NG

    funct var dims: 
    Id          9 x NP
    DV          9 x NP
    Bob_up      9 x NP
%}
% Viscoelastic Constants
    Wi = Qconst.Wi;         % Weissenberg number
    
    d = 3;
    NP = size(Xp,1);
    NG = prod(grid.n+1);    % (nx+1)*(ny+1)*(nz+1)
    w = reshape(w,[d NG]);
    
    Id = repmat([1 0 0 0 1 0 0 0 1]',1,NP);
    Dw = DV_comp(grid,Xp,w);
    Bob_up = dt*( mult_3d(Dw,Bob_k) + multT_3d(Bob_k,Dw) + (1/Wi)*(Id-Bob_k) ) + Bob_k;
end

function [Bup] = B_comp(Fup,Bob_up)
% input: Bob_up, Fup    9 x NP
% output: Bup           9 x NP
    Jup = det_3d(Fup);                      % 1xNP
    detBob = det_3d(Bob_up) + eps;          % 1xNP
    dk = ( (Jup.*Jup) ./ detBob ).^(1/3);   % 1xNP
    Bup = bsxfun(@times,dk,Bob_up);
end

% need: 
function [detA] = det_3d(A)
% A     9 x n
% detA  1 x n
detA = A(1,:).*(A(5,:).*A(9,:)-A(8,:).*A(6,:)) - A(4,:).*(A(2,:).*A(9,:)-A(8,:).*A(3,:)) + A(7,:).*(A(2,:).*A(6,:)-A(5,:).*A(3,:));
end

function [AinvT] = invT_3d(A)
% A         9 x n
% AinvT     9 x n
% detA      1 x n
detA = A(1,:).*(A(5,:).*A(9,:)-A(8,:).*A(6,:)) - A(4,:).*(A(2,:).*A(9,:)-A(8,:).*A(3,:)) + A(7,:).*(A(2,:).*A(6,:)-A(5,:).*A(3,:));
detA = detA + eps;  % avoid dividing by zero

AinvT = [   ( A(5,:).*A(9,:) - A(8,:).*A(6,:) ) ./ detA; 
            ( A(7,:).*A(6,:) - A(4,:).*A(9,:) ) ./ detA; 
            ( A(4,:).*A(8,:) - A(7,:).*A(5,:) ) ./ detA;
            ( A(8,:).*A(3,:) - A(2,:).*A(9,:) ) ./ detA;
            ( A(1,:).*A(9,:) - A(7,:).*A(3,:) ) ./ detA;
            ( A(7,:).*A(2,:) - A(1,:).*A(8,:) ) ./ detA;
            ( A(2,:).*A(6,:) - A(5,:).*A(3,:) ) ./ detA;
            ( A(4,:).*A(3,:) - A(1,:).*A(6,:) ) ./ detA;
            ( A(1,:).*A(5,:) - A(4,:).*A(2,:) ) ./ detA;  ];
end

% maybe need:
function [Ainv] = inv_3d(A)
% A         9 x n
% Ainv      9 x n
% detA      1 x n
detA = A(1,:).*(A(5,:).*A(9,:)-A(8,:).*A(6,:)) - A(4,:).*(A(2,:).*A(9,:)-A(8,:).*A(3,:)) + A(7,:).*(A(2,:).*A(6,:)-A(5,:).*A(3,:));
detA = detA + eps;  % avoid dividing by zero

Ainv = [    ( A(5,:).*A(9,:) - A(8,:).*A(6,:) ) ./ detA;    % 1
            ( A(8,:).*A(3,:) - A(2,:).*A(9,:) ) ./ detA;    % 2
            ( A(2,:).*A(6,:) - A(5,:).*A(3,:) ) ./ detA;    % 3
            ( A(7,:).*A(6,:) - A(4,:).*A(9,:) ) ./ detA;    % 4
            ( A(1,:).*A(9,:) - A(7,:).*A(3,:) ) ./ detA;    % 5
            ( A(4,:).*A(3,:) - A(1,:).*A(6,:) ) ./ detA;    % 6
            ( A(4,:).*A(8,:) - A(7,:).*A(5,:) ) ./ detA;    % 7
            ( A(7,:).*A(2,:) - A(1,:).*A(8,:) ) ./ detA;    % 8
            ( A(1,:).*A(5,:) - A(4,:).*A(2,:) ) ./ detA;  ];% 9
end

function [Atran] = tran_3d(A)
% A         9 x n
% Atran     9 x n

Atran = [A(1,:); A(4,:); A(7,:); A(2,:); A(5,:); A(8,:); A(3,:); A(6,:); A(9,:)];
end

% vectorized (and FASTER than C):
function [M] = mult_3d(A,B)
% M,A,B:     9 x n

% M = A * B
M = [A(1,:).*B(1,:) + A(4,:).*B(2,:) + A(7,:).*B(3,:);
     A(2,:).*B(1,:) + A(5,:).*B(2,:) + A(8,:).*B(3,:);
     A(3,:).*B(1,:) + A(6,:).*B(2,:) + A(9,:).*B(3,:);
     A(1,:).*B(4,:) + A(4,:).*B(5,:) + A(7,:).*B(6,:);
     A(2,:).*B(4,:) + A(5,:).*B(5,:) + A(8,:).*B(6,:);
     A(3,:).*B(4,:) + A(6,:).*B(5,:) + A(9,:).*B(6,:);
     A(1,:).*B(7,:) + A(4,:).*B(8,:) + A(7,:).*B(9,:);
     A(2,:).*B(7,:) + A(5,:).*B(8,:) + A(8,:).*B(9,:);
     A(3,:).*B(7,:) + A(6,:).*B(8,:) + A(9,:).*B(9,:)];
end

function [M] = multT_3d(A,B)
% M,A,B:     9 x n

% M = A * B'
M = [A(1,:).*B(1,:) + A(4,:).*B(4,:) + A(7,:).*B(7,:);
     A(2,:).*B(1,:) + A(5,:).*B(4,:) + A(8,:).*B(7,:);
     A(3,:).*B(1,:) + A(6,:).*B(4,:) + A(9,:).*B(7,:);
     A(1,:).*B(2,:) + A(4,:).*B(5,:) + A(7,:).*B(8,:);
     A(2,:).*B(2,:) + A(5,:).*B(5,:) + A(8,:).*B(8,:);
     A(3,:).*B(2,:) + A(6,:).*B(5,:) + A(9,:).*B(8,:);
     A(1,:).*B(3,:) + A(4,:).*B(6,:) + A(7,:).*B(9,:);
     A(2,:).*B(3,:) + A(5,:).*B(6,:) + A(8,:).*B(9,:);
     A(3,:).*B(3,:) + A(6,:).*B(6,:) + A(9,:).*B(9,:)];
end

function [M] = Tmult_3d(A,B)
% M,A,B:     9 x n

% M = A' * B
M = [A(1,:).*B(1,:) + A(2,:).*B(2,:) + A(3,:).*B(3,:);
     A(4,:).*B(1,:) + A(5,:).*B(2,:) + A(6,:).*B(3,:);
     A(7,:).*B(1,:) + A(8,:).*B(2,:) + A(9,:).*B(3,:);
     A(1,:).*B(4,:) + A(2,:).*B(5,:) + A(3,:).*B(6,:);
     A(4,:).*B(4,:) + A(5,:).*B(5,:) + A(6,:).*B(6,:);
     A(7,:).*B(4,:) + A(8,:).*B(5,:) + A(9,:).*B(6,:);
     A(1,:).*B(7,:) + A(2,:).*B(8,:) + A(3,:).*B(9,:);
     A(4,:).*B(7,:) + A(5,:).*B(8,:) + A(6,:).*B(9,:);
     A(7,:).*B(7,:) + A(8,:).*B(8,:) + A(9,:).*B(9,:)];
end

% ======================== Wrapper Functions =============================
function fp = G2P_PIC(grid,Xp,fg,b)
%{  
var dimensions:
    Xp:     NP x 3
    fg:     1 x NG
    b =     [bx,by,bz]

    fp:     1 x NP  
%}

    % Parse inputs (unpack)
    Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3); 
    hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
    nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
    Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1); 
    bx = b(1);      by = b(2);      bz = b(3);

    % Call C function
    fp = G2P_PIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, fg, bx,by,bz);
end

function fg = P2G_PIC(grid,Xp,fp,b)
%{  
var dimensions:
    Xp:     NP x 3
    fp:     1 x NP
    b =     [bx,by,bz]

    fg:     1 x NG  
%}

    % Parse inputs (unpack)
    Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3); 
    hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
    nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
    Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1); 
    bx = b(1);      by = b(2);      bz = b(3);

    % Call C function
    fg = P2G_PIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, fp, bx,by,bz);
end

function [vp] = G2P_vec(grid,Xp,vg,b)
%{  
var dimensions:
    Xp:     NP x 3
    vg:     3 x NG
    b =     [bx,by,bz]

    vp:     3 x NP  
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1);
fg1 = vg(1,:);  fg2 = vg(2,:);  fg3 = vg(3,:);
bx = b(1);      by = b(2);      bz = b(3);

% Call C function
[fp1,fp2,fp3] = G2P_vec_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, fg1,fg2,fg3, bx,by,bz);

% Package Outputs
vp = [fp1; fp2; fp3];
end

function [vg] = P2G_vec(grid,Xp,vp,b)
%{  
var dimensions:
    Xp:     NP x 3
    vp:     3 x NP
    b =     [bx,by,bz]

    vg:     3 x NG  
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1);
fp1 = vp(1,:);  fp2 = vp(2,:);  fp3 = vp(3,:);
bx = b(1);      by = b(2);      bz = b(3);

% Call C function
[fg1,fg2,fg3] = P2G_vec_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, fp1,fp2,fp3, bx,by,bz);

% Package Outputs
vg = [fg1; fg2; fg3];
end

function [Fup] = F_comp(grid,Xp,dt,v,F)
%{  
var dimensions:
    Xp:     NP x 3
    dt:     scalar
    v:      3 x NG
    F:      9 x NP
    Fi:     1 x NP

    Fup:    9 x NP
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1);
vx = v(1,:);    vy = v(2,:);    vz = v(3,:);
F1 = F(1,:);    F2 = F(2,:);    F3 = F(3,:);
F4 = F(4,:);    F5 = F(5,:);    F6 = F(6,:);
F7 = F(7,:);    F8 = F(8,:);    F9 = F(9,:);

% Call C function
[Fup1,Fup2,Fup3,Fup4,Fup5,Fup6,Fup7,Fup8,Fup9] = F_comp_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, dt, vx,vy,vz, F1,F2,F3,F4,F5,F6,F7,F8,F9);

% Package Outputs
Fup = [Fup1; Fup2; Fup3; Fup4; Fup5; Fup6; Fup7; Fup8; Fup9];
end

function [Phi,Vp,B] = G2P_APIC(grid,Xp,xup,w)
%{  
var dimensions:
    Xp:     NP x 3
    xup,w:  3 x NG
    
    Phi:    NP x 3
    Vp:     3 x NP
    B:      9 x NP
%}

% Parse inputs (unpack)
Ox = grid.O(1);     Oy = grid.O(2);     Oz = grid.O(3);
hx = grid.h(1);     hy = grid.h(2);     hz = grid.h(3);
nx = grid.n(1);     ny = grid.n(2);     nz = grid.n(3);
Yp = Xp(:,2);       Zp = Xp(:,3);       Xp = Xp(:,1);
xup1 = xup(1,:);    xup2 = xup(2,:);    xup3 = xup(3,:);
w1 = w(1,:);        w2 = w(2,:);        w3 = w(3,:);

% Call C function
[Phi1,Phi2,Phi3, Vp1,Vp2,Vp3, B1,B2,B3,B4,B5,B6,B7,B8,B9] = G2P_APIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, xup1,xup2,xup3, w1,w2,w3);

% Package Outputs
Phi = [Phi1, Phi2, Phi3];
Vp = [Vp1; Vp2; Vp3];
B = [B1; B2; B3; B4; B5; B6; B7; B8; B9;];
end

function [mom_g,mass_g] = P2G_APIC(grid,Xp,B,mom_p,mass_p)
%{  
var dimensions:
    Xp:         NP x 3
    B:          9 x NP
    mom_p:      3 x NP
    mass_p:     NP x 1
    
    mom_g:      3 x NG
    mass_g:     NG x 1
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1);
B1 = B(1,:);    B2 = B(2,:);    B3 = B(3,:); 
B4 = B(4,:);    B5 = B(5,:);    B6 = B(6,:);
B7 = B(7,:);    B8 = B(8,:);    B9 = B(9,:); 
mom_p1 = mom_p(1,:); mom_p2 = mom_p(2,:); mom_p3 = mom_p(3,:);
mass_p = mass_p(:);

% Call C function
[mom_g1,mom_g2,mom_g3, mass_g] = P2G_APIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, B1,B2,B3,B4,B5,B6,B7,B8,B9, mom_p1,mom_p2,mom_p3, mass_p);

% Package Outputs
mom_g = [mom_g1; mom_g2; mom_g3];
end

function [DV] = DV_comp(grid,Xp,v)
%{  
var dimensions:
    Xp:     NP x 3
    v:      3 x NG

    DV:    9 x NP
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2);   Zp = Xp(:,3);   Xp = Xp(:,1);
vx = v(1,:);    vy = v(2,:);    vz = v(3,:);

% Call C function
[DV1,DV2,DV3,DV4,DV5,DV6,DV7,DV8,DV9] = DV_comp_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, vx,vy,vz);

% Package Outputs
DV = [DV1; DV2; DV3; DV4; DV5; DV6; DV7; DV8; DV9];
end

