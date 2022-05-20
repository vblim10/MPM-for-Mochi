% mochi_slice.m 
% Slicing Mochi Ice Cream
% Pang-thao (butcher knife) & "cutting board"

% *** IMPORT WORKSPACE *** wkspc_mochi_ctr_xy_8_8.mat ***

d = 3;              % dimension
CFL = 0.25; 
tf  = 1.0;
g = 1.0;            % gravitational constant
dtk = 0.0166;       % (dtk=0.0333, 30 fps) (dtk=0.0667, 15 fps)

% Material & Plot settings
ppc  = 2;           % # particles per grid cell = ppc^d
psize = 25;         % for plotting particles
sqSz = 100;         % for plotting grid nodes

% ====================== VISCOELASTIC CONSTANTS ==========================
% Ice Cream
md_p_cream = 1;             % particle mass density
mu_cream = 1.5;             % shear modulus
La_cream = 15;              % Lame's 1st parameter
mu_N_cream = 1.0*10^-2;     % viscous Newton constant
Wi_cream = 0.30;            % Weissenberg number

% Mochi
md_p_mochi = 1.28;          % particle mass density
mu_mochi = 6;               % shear modulus
La_mochi = 60;              % Lame's 1st parameter
mu_N_mochi = 1.0*10^-2;     % viscous Newton constant
Wi_mochi = 3;               % Weissenberg number

% ======================= LAGRANGIAN PARTICLES ===========================
% *** IMPORT WORKSPACE *** (Phi_up, Idx_cream, Idx_mochi)
% Particle postions for Ice cream & Mochi
Phi_k = Phi_up;   
NP = size(Phi_k,1);             % total # particles
NP_cream = length(Idx_cream);   % # ice cream particles
NP_mochi = length(Idx_mochi);   % # mochi particles

% Initial velocity
v1x = 0.0;  v1y = 0.0; v1z = 0.0;   % vel = <0,0,0> (mochi ice cream)

% VE params struct for funct. handles [ice cream, mochi]
Qconst.mu = [mu_cream*ones(1,NP_cream),     mu_mochi*ones(1,NP_mochi)];          
Qconst.La = [La_cream*ones(1,NP_cream),     La_mochi*ones(1,NP_mochi)];          
Qconst.mu_N = [mu_N_cream*ones(1,NP_cream), mu_N_mochi*ones(1,NP_mochi)];      
Qconst.Wi = [Wi_cream*ones(1,NP_cream),     Wi_mochi*ones(1,NP_mochi)];  

% ======================== EULERIAN MESH ================================= 
h = [1.0, 1.0, 1.0];            % grid stepsize  
% h = [0.5, 0.5, 0.5];            % finer mesh

% grid struct
grid.h = h;
[grid,Xg,NG] = EulGrid(Phi_k,grid); % NG = # grid nodes
                                    % Xg = grid node positions (NGx3)
    
% ======================== KINEMATIC BLADE =============================== 
% Blade constants
hb = 0.32;  % blade node spacing & thickness
a = 31;     % ellipse long "radius" (scaled with hb)

% Default Blade
[blade,Phi_blade,n_co] = BladeLevelSet(a,hb,grid);

% Edit Blade
blade.d = [0 -1];           % direction of blade (in xz-plane)
blade.o = [8.0 23.0];       % x,z-coord of ellipse center 
blade.x_co = Phi_blade + [blade.o(1) 0 blade.o(2)];     % (NObj x 3)
blade.v_co = [0; 0; -15];   % (3 x 1)   init. vel = <0,0,-15>
blade.n_co = n_co;          % (3 x NObj)
blade.a_co = [0; 0; -0.25]; % (3 x 1) 
    
% ========================== CUTTING BOARD =============================== 
[Xi,Xj,Xk] = meshgrid((1:h(1):15),(1:h(2):15),-0.1);
NP_board = length(Xi(:));
board.x_co = [Xi(:),Xj(:),Xk(:)];           % (NObj x 3)
board.v_co = [0; 0; 0];                     % (3 x 1)   velocity = <0,0,0>     
board.n_co = repmat([0 0 1]',1,NP_board);   % (3 x NObj)

% =========================== CHECK GRAPH ================================
figure(1)
scatter3(Phi_k(Idx_cream,1),Phi_k(Idx_cream,2),Phi_k(Idx_cream,3),psize,'filled','m'); hold on; % ice cream
scatter3(Phi_k(Idx_mochi,1),Phi_k(Idx_mochi,2),Phi_k(Idx_mochi,3),psize,'filled','g'); hold on; % mochi
scatter3(blade.x_co(:,1),blade.x_co(:,2),blade.x_co(:,3),psize,'filled','k'); hold on; % blade
scatter3(board.x_co(:,1),board.x_co(:,2),board.x_co(:,3),psize,'filled','k'); hold on; % board
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','b'); hold on; % grid
hold off;

% ============================ Step k = 0 ================================
tic
k = 0;      
t = 0;      % time at step k

% Compute Particle Volumes (Tot_Volume/NP)
    Vol = prod(h)*(1/ppc^d)*ones(NP,1); 
    
% Set Particle Mass (NPx1)
    mass_p = [ md_p_cream*ones(NP_cream,1); 
               md_p_mochi*ones(NP_mochi,1)].*Vol;  
    
% MPM Initial Values (k=0)
    % Deformation gradient (F = Id for all particles) 
    Id = eye(d);
    Fk = repmat(Id(:),[1,NP]);      % 9xNP
    Jk = det_3d(Fk);                % 1xNP
    
    % Elastic Stress vars
    Bob_k = repmat(Id(:),[1,NP]);   % 9xNP
    Bk = B_comp(Fk,Bob_k);          % 9xNP
    
    % APIC var
    Bm = Fk;                        % 9xNP
    
    % Particle Velocity (3xNP)
    Uk = [v1x*ones(1,NP); 
          v1y*ones(1,NP);
          v1z*ones(1,NP)];     
       
    % Tolerance (for minimization)
    tol = min([max(h)^4, 1e-4, max(h)^-4]); % cases: h<1, h=1, h>1
    
    % "Actual" Energy over time 
    % funct. handle: grid, Vol, mass_p, Qconst, & g are constant over time
    RealEngy = @(J,B_El,Uk,w,Xp) RealEngy_comp(J,B_El,Uk,w,Xp,grid,Vol,mass_p,Qconst,g);
    u = P2G_vec(grid,Phi_k,Uk,[0 0 0]); w = u; % so we can compute init. Engy
    
    % ANIMATION
    T = [t,t]; 
    PHX = [Phi_k(:,1), Phi_k(:,1)]; 
    PHY = [Phi_k(:,2), Phi_k(:,2)];
    PHZ = [Phi_k(:,3), Phi_k(:,3)];
    JPK = [det_3d(Fk)', det_3d(Fk)'];
    BX = [blade.x_co(:,1), blade.x_co(:,1)];
    BY = [blade.x_co(:,2), blade.x_co(:,2)];
    BZ = [blade.x_co(:,3), blade.x_co(:,3)];

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
    w = reshape(w,[d NG]);  % 3xNG
    
    % COLLISION HANDLING: grid velocities
    [Idx1] = CollisionDetect_Blade(Xg,blade);
    [Idx2] = CollisionDetect_Board(Xg,board);
    [v_up1] = CollisionHandling(Idx1,Xg(Idx1,:),w(:,Idx1),blade);
    [v_up2] = CollisionHandling(Idx2,Xg(Idx2,:),w(:,Idx2),board);
    w(:,Idx1) = v_up1; w(:,Idx2) = v_up2;
    
    % APIC TRANSFER grid to particle (advection) 
    xup = Xg' + dtk*w;
    [Phi_up,Vk,Bm] = G2P_APIC(grid,Phi_k,xup,w); 
    
    % COLLISION HANDLING: particle velocities
    [Idx1] = CollisionDetect_Blade(Phi_up,blade); 
    [Idx2] = CollisionDetect_Board(Phi_up,board);
    [v_up1] = CollisionHandling(Idx1,Phi_up(Idx1,:),Vk(:,Idx1),blade);
    [v_up2] = CollisionHandling(Idx2,Phi_up(Idx2,:),Vk(:,Idx2),board);
    Vk(:,Idx1) = v_up1; Vk(:,Idx2) = v_up2;
    
    % Particle Advection (collision corrected Lag velocity)
    Phi_up = Phi_k + dtk*Vk';

    % Graph 
    scatter3(Phi_up(Idx_cream,1),Phi_up(Idx_cream,2),Phi_up(Idx_cream,3),psize,'filled','m'); hold on; % ice cream
    scatter3(Phi_up(Idx_mochi,1),Phi_up(Idx_mochi,2),Phi_up(Idx_mochi,3),psize,'filled','g'); hold on; % mochi
    scatter3(blade.x_co(:,1),blade.x_co(:,2),blade.x_co(:,3),psize,'filled','k'); hold on; % blade
    scatter3(board.x_co(:,1),board.x_co(:,2),board.x_co(:,3),psize,'filled','k'); hold on; % board
    scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on; % grid
    title(num2str(kb),num2str(t)); xlabel(num2str(k)); ylabel(num2str(dtk));
    xlim([0 16]); ylim([0 16]); zlim([-5 35]);
    az = 0.0447; el = 0.0677;   % side view (xz-plane)
    view([az,el]);  
    hold off; pause(1e-10);
    
    % ANIMATION
    T = [T,t]; 
    PHX = [PHX, Phi_up(:,1)]; 
    PHY = [PHY, Phi_up(:,2)]; 
    PHZ = [PHZ, Phi_up(:,3)];
    JPK = [JPK, det_3d(Fk)'];
    BX = [BX, blade.x_co(:,1)];
    BY = [BY, blade.x_co(:,2)];
    BZ = [BZ, blade.x_co(:,3)];
    
    % Update vars for next time step 
    k = k + 1;
    t = t + dtk;
    Fk = F_comp(grid,Phi_k,dtk,w,Fk);   % 9xNP
    Jk = det_3d(Fk);                    % 1xNP
    Bob_k = Bob_comp(Bob_k,dtk,Qconst,grid,Phi_k,w); % 9xNP
    Bk = B_comp(Fk,Bob_k);              % 9xNP
    Phi_k = Phi_up;                     % NPx3
    Uk = Vk;                            % 3xNP
    
    % Update Blade 
    [blade] = BladeUp(blade,dtk);
    
    % Update Eulerian grid
    [grid,Xg,NG] = EulGrid(Phi_k,grid);
    
end
time = toc;

% "Actual" Energy over time 
EvT(k+1,:) = [t,RealEngy(Jk,Bk,Uk,w,Phi_k)]; % [t,E,Q,K,Q_El,Q_N]

% ++++++++++++++++++++++++++++ PLOTS +++++++++++++++++++++++++++++++++++++

% Energy Plots
figure(2)
plot(EvT(:,1),EvT(:,2),'r');   hold on; % Total Energy
plot(EvT(:,1),EvT(:,3),'-go'); hold on; % Potential Energy
plot(EvT(:,1),EvT(:,4),'-bs'); hold on; % Kinetic Energy
legend('Total Energy','Potential Energy','Kinetic Energy');
xlabel('time t'); title('Energy'); hold off;

figure(3)
plot(EvT(:,1),EvT(:,2),'r');    % Total Energy
xlabel('time t'); title('Total Energy'); 

% ------------------------- ANIMATION DATA -------------------------------
% Write Energy data to text file: EvT = [t,E,Q,K,Q_El,Q_N,G]
writematrix(EvT,'Engy.dat','Delimiter',';');  % (k+1 x 7)

% Write Particle positions to text file
writematrix(PHX,'PHX.dat','Delimiter',';');  % (x-coords) (NP x k+1)
writematrix(PHY,'PHY.dat','Delimiter',';');  % (y-coords) (NP x k+1)
writematrix(PHZ,'PHZ.dat','Delimiter',';');  % (z-coords) (NP x k+1)

% Write Blade (Collider) node positions to text file 
writematrix(BX,'BX.dat','Delimiter',';');  % (x-coords) (NObj x k+1)
writematrix(BY,'BY.dat','Delimiter',';');  % (y-coords) (NObj x k+1)
writematrix(BZ,'BZ.dat','Delimiter',';');  % (z-coords) (NObj x k+1)

% Write J=det(Fp) particle data to text file
writematrix(JPK,'Jpk.dat','Delimiter',';');  % (NP x k+1)

% -------------------------- ANIMATED SIM --------------------------------
% SAVE Workspace and IMPORT into mochi_slice_ani.m  
% in order to render simulation as gif/mp4

% +++++++++++++++++++++++++++ FUNCTIONS ++++++++++++++++++++++++++++++++++
function [grid,Xg,NG] = EulGrid(Phi_k,grid)
% function to construct Eulerian Grid (rect. box) & update grid struct

    % unpack grid struct
    h = grid.h;     % 1x3

    % Find range of Lagrangian material 
    xmin = round(min(Phi_k(:,1))); xmax = round(max(Phi_k(:,1)));
    ymin = round(min(Phi_k(:,2))); ymax = round(max(Phi_k(:,2)));
    zmin = round(min(Phi_k(:,3))); zmax = round(max(Phi_k(:,3)));
    
    % Construct Eulerian mesh (+2h padding)
    x = xmin-2*h(1):h(1):xmax+2*h(1);   
    y = ymin-2*h(2):h(2):ymax+2*h(2);
    z = zmin-2*h(3):h(3):zmax+2*h(3);
    nx = length(x)-1;                   % # x-axis partitions
    ny = length(y)-1;                   % # y-axis partitions
    nz = length(z)-1;                   % # z-axis partitions
    [Ei,Ej,Ek] = meshgrid(x,y,z);       % grid node positions
    
    NG = length(Ei(:));             % # grid nodes
    Xg = [Ei(:),Ej(:),Ek(:)];       % grid node positions (NG x 3)
    O = [min(x); min(y); min(z)];   % "origin" = bottom-left corner of grid
    
    % grid struct
    grid.O = O';                        % 1x3
    grid.h = h;                         % 1x3
    grid.n = [nx, ny, nz];              % 1x3
    grid.height = Ek(:);                % NGx1
    grid.x = x; grid.y = y; grid.z = z; % row vectors
end

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

    funct var dims: 
    Vol_k               NP x 1
    trB_El              1 x NP
    H                   1 x NP
    DV,S                9 x NP
    E,Q,K,Q_El,Q_N      scalar
    Engy                1 x 6
%}

% dimension
    d = 3;

% Elastic & Viscous constants (1xNP)
    mu = Qconst.mu;     % shear modulus
    La = Qconst.La;     % Lame's 1st parameter
    mu_N = Qconst.mu_N; % viscous Newton constant
    
% Potential Energy 
    Vol_k = (J').*Vol;
    trB = Bk(1,:) + Bk(5,:) + Bk(9,:);
    % H(B_El)
    H = 0.5*mu.*(trB - d) - mu.*log(J) + 0.5*La.*(J-1).*(J-1);
    Q_El = dot(H,Vol');     % elastoplastic portion of potential energy
    DV = DV_comp(grid,Xp,w);
    S = 0.5*(DV + tran_3d(DV));
    Q_N = ( mu_N.* dot(S,S,1) )*Vol_k;  % Newtonian viscous potential energy
    G = g*dot(mass_p, Xp(:,3));         % Gravitational potential energy 
    Q = Q_El + Q_N + G;                 % Potential Energy
    
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

% Elastic & Viscous constants (1xNP)
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
    Wi_inv = bsxfun(@rdivide,1,Wi); % 1/Wi (1xNP)
    Bob_up = Bob_k + dt*( mult_3d(Dw,Bob_k) + multT_3d(Bob_k,Dw) );
    Bob_up = Bob_up + dt*bsxfun(@times,Wi_inv,Id-Bob_k);
    detBob = det_3d(Bob_up) + eps;          % 1xNP
    dk = ( (Jup.*Jup) ./ detBob ).^(1/3);   % 1xNP
    Bup = bsxfun(@times,dk,Bob_up);         
    
% Elastic Energy: El
    trB = Bup(1,:) + Bup(5,:) + Bup(9,:);                   % 1xNP
    PB = 0.5*mu.*(trB - d) - mu.*log(Jup) + La.*(Jup-1).^2; % 1xNP
    El = dot(PB,Vol');   
    
% Viscous Energy: En
    Dw_sym = 0.5*( Dw + tran_3d(Dw) );  % 9xNP
    R = sum(Dw_sym.*Dw_sym,1);          % 1xNP
    Vol_up = (Jup').*Vol;               % NPx1
    En = dot( mu_N.*R ,Vol_up');

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
    dk          1 x NP
    detBob      1 x NP
    Bup         9 x NP
    P,P1,P2     9 x NP
    trP         1 x NP
    Bob_up_invT 9 x NP
    Dw_sym      9 x NP
    Q           9 x NP
    gradPot, gradK, gradE   (3*NG) x 1
%}

% Elastic & Viscous constants (1xNP)
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
    Wi_inv = bsxfun(@rdivide,1,Wi); % 1/Wi (1xNP)
    Bob_up = Bob_k + dt*( mult_3d(Dw,Bob_k) + multT_3d(Bob_k,Dw) );
    Bob_up = Bob_up + dt*bsxfun(@times,Wi_inv,Id-Bob_k);
    detBob = det_3d(Bob_up) + eps;              % 1xNP
    dk = ( (Jup.*Jup) ./ detBob ).^(1/3);       % 1xNP
    Bup = bsxfun(@times,dk,Bob_up);             
    
% Elastic Energy: gradEl
    P = 0.5*bsxfun(@times,mu,Bup); 
    P = P + 0.5*bsxfun(@times,( La.*(Jup.*Jup - Jup)- mu ),Id); % 9xNP
    % gradE_OB_1
    trP = P(1,:) + P(5,:) + P(9,:);             % 1xNP
    Bob_up_invT = invT_3d(Bob_up);              % 9xNP
    P1 = 2*multT_3d(invT_3d(Fup),Fk) - multT_3d(Bob_up_invT,Bob_k) - Tmult_3d(Bob_k,Bob_up_invT);
    P1 = (dt/3) * bsxfun(@times,trP,P1);        % 9xNP
    % gradE_OB_2
    P2 = 2*dt*P; 
    
% Viscous Energy: gradEn
    Dw_sym = 0.5*( Dw + tran_3d(Dw) );          % 9xNP
    Q = 2*bsxfun(@times,mu_N,Dw_sym);           % 9xNP
    
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
% Viscoelastic Constants (1xNP)
    Wi = Qconst.Wi;         % Weissenberg number
    
    d = 3;
    NP = size(Xp,1);
    NG = prod(grid.n+1);    % (nx+1)*(ny+1)*(nz+1)
    w = reshape(w,[d NG]);
    
    Id = repmat([1 0 0 0 1 0 0 0 1]',1,NP);
    Dw = DV_comp(grid,Xp,w);
    
    Wi_inv = bsxfun(@rdivide,1,Wi); % 1/Wi (1xNP)
    Bob_up = Bob_k + dt*( mult_3d(Dw,Bob_k) + multT_3d(Bob_k,Dw) );
    Bob_up = Bob_up + dt*bsxfun(@times,Wi_inv,Id-Bob_k);
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

% ========================= Wrapper Functions ============================
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
    Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1); 
    bx = b(1); by = b(2); bz = b(3);

    % Call C function
    fp = G2P_PIC_3d(Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, Xp, Yp, Zp, fg, bx, by, bz);
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
    Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1); 
    bx = b(1); by = b(2); bz = b(3);

    % Call C function
    fg = P2G_PIC_3d(Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, Xp, Yp, Zp, fp, bx, by, bz);
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
Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1);
fg1 = vg(1,:); fg2 = vg(2,:); fg3 = vg(3,:);
bx = b(1); by = b(2); bz = b(3);

% Call C function
[fp1, fp2, fp3] = G2P_vec_3d(Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, Xp, Yp, Zp, fg1, fg2, fg3, bx, by, bz);

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
Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1);
fp1 = vp(1,:); fp2 = vp(2,:); fp3 = vp(3,:);
bx = b(1); by = b(2); bz = b(3);

% Call C function
[fg1, fg2, fg3] = P2G_vec_3d(Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, Xp, Yp, Zp, fp1, fp2, fp3, bx, by, bz);

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
Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1);
% dt = dt;
vx = v(1,:); vy = v(2,:); vz = v(3,:);
F1 = F(1,:); F2 = F(2,:); F3 = F(3,:);
F4 = F(4,:); F5 = F(5,:); F6 = F(6,:);
F7 = F(7,:); F8 = F(8,:); F9 = F(9,:);

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
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1);
xup1 = xup(1,:); xup2 = xup(2,:); xup3 = xup(3,:);
w1 = w(1,:); w2 = w(2,:); w3 = w(3,:);

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
Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1);
B1 = B(1,:); B2 = B(2,:); B3 = B(3,:); 
B4 = B(4,:); B5 = B(5,:); B6 = B(6,:);
B7 = B(7,:); B8 = B(8,:); B9 = B(9,:); 
mom_p1 = mom_p(1,:); mom_p2 = mom_p(2,:); mom_p3 = mom_p(3,:);
mass_p = mass_p(:);

% Call C function
[mom_g1,mom_g2,mom_g3,mass_g] = P2G_APIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, B1,B2,B3,B4,B5,B6,B7,B8,B9, mom_p1,mom_p2,mom_p3, mass_p);

% Package Outputs
mom_g = [mom_g1; mom_g2; mom_g3];
end

function [DV] = DV_comp(grid,Xp,v)
%{  
var dimensions:
    Xp:     NP x 3
    v:      3 x NG

    DV:     9 x NP
%}

% Parse inputs (unpack)
Ox = grid.O(1); Oy = grid.O(2); Oz = grid.O(3);
hx = grid.h(1); hy = grid.h(2); hz = grid.h(3);
nx = grid.n(1); ny = grid.n(2); nz = grid.n(3);
Yp = Xp(:,2); Zp = Xp(:,3); Xp = Xp(:,1);
vx = v(1,:); vy = v(2,:); vz = v(3,:);

% Call C function
[DV1,DV2,DV3,DV4,DV5,DV6,DV7,DV8,DV9] = DV_comp_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, vx,vy,vz);

% Package Outputs
DV = [DV1; DV2; DV3; DV4; DV5; DV6; DV7; DV8; DV9];
end

% =========================== RIGID OBJECTS ==============================
function [blade,Phi_blade,n_co] = BladeLevelSet(a,hb,grid)
%{
function to compute (default) Level Set nodes and normals for the blade 
center:(0,0,0), front & back: yz-plane, vertical in xz-plane

    NObj = # object nodes

Input:  a           long "radius" of blade (scaled with hb)
        hb          node spacing & blade thickness 
        grid        (struct) Eul. grid

Output: blade       (struct)  
        Phi_blade   node positions
        n_co        node normals
%}

% Default Blade
blade.d = [0 -1];                   % direction of blade1 (in xz-plane)
blade.o = [0 0];                    % x,z-coord of ellipse center 
blade.a = a*hb;                     % ellipse long "radius" 
blade.b = hb/2;                     % ellipse short "radius"
blade.y = [min(grid.y)+0.5*hb: hb: max(grid.y)-0.5*hb]'; % node y-coords
blade.ny = length(blade.y);         % # nodes along y-axis

xmin = blade.o(1) - blade.b;
xmax = blade.o(1) + blade.b;
zmin = blade.o(2) - blade.a;
zmax = blade.o(2) + blade.a;
    
% Blade node positions
%                          x            y        z
    [Xi,Xj,Xk] = meshgrid([xmin, xmax], blade.y, [zmin + hb: hb: zmax - hb]);
    % Indices for orange nodes (meshgrid orders nodes: y,x,z)
    orgl = zeros(size(Xi,1),size(Xi,2),size(Xi,3));
    orgl(1,:,:)=1; orgl = find(orgl(:));                % left
    orgr = zeros(size(Xi,1),size(Xi,2),size(Xi,3));
    orgr(end,:,:)=1; orgr = find(orgr(:));              % right
    % Indices for cyan nodes
    cynf = zeros(size(Xi,1),size(Xi,2),size(Xi,3));
    cynf(2:end-1,1,:) = 1; cynf = find(cynf(:));        % front
    cynb = zeros(size(Xi,1),size(Xi,2),size(Xi,3));
    cynb(2:end-1,end,:) = 1; cynb = find(cynb(:));      % back
    % Initial Positions
    Phi_blade = [Xi(:),                       Xj(:),   Xk(:);                  % mid
                 blade.o(1)*ones(blade.ny,1), blade.y, zmax*ones(blade.ny,1);  % top 
                 blade.o(1)*ones(blade.ny,1), blade.y, zmin*ones(blade.ny,1)]; % bot 
    NP_blade = size(Phi_blade,1);
    % Indices for green nodes
    grnb = NP_blade-blade.ny+1: 1: NP_blade;        % bot
    grnt = NP_blade-(2*blade.ny)+1: 1: min(grnb)-1; % top
    
    % Check indices
figure(3)
scatter3(Phi_blade(orgl,1),Phi_blade(orgl,2),Phi_blade(orgl,3),'filled','y'); hold on;  % orgl
scatter3(Phi_blade(orgr,1),Phi_blade(orgr,2),Phi_blade(orgr,3),'filled','m'); hold on;  % orgr
scatter3(Phi_blade(grnt,1),Phi_blade(grnt,2),Phi_blade(grnt,3),'filled','g'); hold on;  % grnt
scatter3(Phi_blade(grnb,1),Phi_blade(grnb,2),Phi_blade(grnb,3),'filled','r'); hold on;  % grnb
scatter3(Phi_blade(cynf,1),Phi_blade(cynf,2),Phi_blade(cynf,3),'filled','c'); hold on;  % cynf
scatter3(Phi_blade(cynb,1),Phi_blade(cynb,2),Phi_blade(cynb,3),'filled','b'); hold on;  % cynb
title('Blade Level Set');
hold off;
    
% Blade normals (3 x NP_blade)
n_co = zeros(3,NP_blade);
n_co(:,orgl) = repmat([ 0 -1  0]',1,length(orgl)); % orgl
n_co(:,orgr) = repmat([ 0  1  0]',1,length(orgr)); % orgr
n_co(:,grnt) = repmat([ 0  0  1]',1,length(grnt)); % grnt
n_co(:,grnb) = repmat([ 0  0 -1]',1,length(grnb)); % grnb
n_co(:,cynf) = repmat([-1  0  0]',1,length(cynf)); % cynf
n_co(:,cynb) = repmat([ 1  0  0]',1,length(cynb)); % cynb
    
% Set blade struct values
blade.x_co = Phi_blade;
blade.n_co = n_co;
end

function [blade] = BladeUp(blade,dt)
% function to update Blade level set (struct)
% *** simple translation (no rotation) ***

% unpack struct
x_co = blade.x_co;
v_co = blade.v_co;
n_co = blade.n_co;
a_co = blade.a_co;

% Update blade node positions
x_co = x_co + dt*v_co'; 

blade.d;  % direction of blade (in xz-plane), stays same (no rotation)

% Update blade center
blade.o = blade.o + dt*[v_co(1) v_co(3)]; % x,z-coord of ellipse ctr (1x2)

% Update node velocities
v_co = v_co + dt*a_co;

blade.v_co = v_co;
blade.x_co = x_co;
blade.n_co = n_co;   % stay same (no rotation)
end

function [Idx] = CollisionDetect_Blade(x,blade)
%{ 
% Collision Detection function of lenticular blade (cross section=ellipse)
  Movement only in the xz-plane, so y-coords are const

Input:  x   (NPoG x 3)  particle/grid node positions
        d   (1 x 2)     direction of blade (ONLY in xz-plane)
        o   (1 x 2)     x,z-coord of center of ellipse
        a,b scalar      long, short widths of ellipse 
        y_min, x_min    domain of blade in y-axis 

Output: Idx (NC x 1)    indices of particle/grid nodes "inside" blade
%}
% Unpack Struct 
    d = blade.d; o = blade.o; a = blade.a; b = blade.b;
    y_min = min(blade.y);
    y_max = max(blade.y);

% Compute: (xz-coords ONLY)
    c = sqrt(a^2 - b^2);        % focal length
    foc1 = o - (c/norm(d,2))*d; % focal point 1     (1 x 2)
    foc2 = o + (c/norm(d,2))*d; % focal point 2     (1 x 2)
    m = bsxfun( @minus, foc1, [x(:,1), x(:,3)] ); % (NPoG x 2)  
    n = bsxfun( @minus, foc2, [x(:,1), x(:,3)] ); % (NPoG x 2)  

% Check if inside blade 
    valid = ( vecnorm(m,2,2) + vecnorm(n,2,2) <= 2*a );        % (NPoG x 1)
    valid = valid .* ( x(:,2) > y_min ) .* ( x(:,2) < y_max ); % (NPoG x 1)
    Idx = find(valid);                                         % (NC x 1)
end

function [Idx] = CollisionDetect_Board(x,board)
%{
Collision Detection function of Cutting Board (static)

Input:  x       (NPoG x 3) particle/grid node positions
        board   (struct)

board:  x_co        (NObj x 3)  node positions
        v_co        (3 x 1)     node velocities
        n_co        (3 x NObj)  normal vectors
        
Output: Idx (NC x 1)    indices of particle/grid nodes "inside" board
%}
% Unpack board struct
    x_co = board.x_co;
    top = max(x_co(:,3));   % z-coord of cutting board surface

% Check if inside cutting board   
    in = (x(:,3) <= top);   % (NPoG x 1) (bool)
    Idx = find(in);         % (NC x 1)
end

% ======================== COLLISION FUNCTIONS ===========================
function [v_up] = CollisionHandling(Idx,x,v,object)
%{
% function to update the velocities of particles or grid nodes 
  "colliding" with the object (that undergoes "rigid" transformations)

  NC = # colliding particles/grid nodes
  NObj = # object nodes

Input:  x       (NC x 3)    particle/grid node positions
        v       (3 x NC)    particle/grid node velocities
        v_co    (3 x 1)     object velocity (same for all nodes)
        x_co    (NObj x 3)  object node positions
        n       (3 x NObj)  object node normals (*ALL nodes are bdry nodes)

Vars:   v_rel   (3 x NC)  relative velocity
        v_n     (1 x NC)  dot(v_rel,n_co)
        v_t     (3 x NC)  tangential component of relative velocity

Output: v_up    (3 x NC)
%}
% if not colliding, return unchanged velocity v
    if(isempty(Idx))
       v_up = v;
       return; 
    end
    
% friction coeff
    mu = 0.42;  % sliding friction coeff. (steel)
%     mu = 0.65;   % static friction coeff. (steel)

% Unpack object struct
    v_co = object.v_co;
    x_co = object.x_co;
    n_co = object.n_co;

% Compute: 
    v_rel = bsxfun(@minus, v, v_co);    % (3 x NC)
    
% For each particle/grid node
    % Find index of closest object node for each particle/grid node
    [Idx_co] = ClosestObjNode(x,x_co);  % (NC x 1)
    
    % Compute:
    v_n = dot(v_rel,  n_co(:,Idx_co), 1);              % (1 x NC)
    v_t = v_rel - bsxfun(@times, n_co(:,Idx_co), v_n); % (3 x NC)
    v_t_norm = vecnorm(v_t,2,1) + eps;                 % (1 x NC)
    
    
    % Update: (collision applied)
    v_rel = v_t + mu * bsxfun(@times, v_t, (v_n ./ v_t_norm) ); % (3 x NC)
    v_up = bsxfun(@plus,v_rel, v_co);                           % (3 x NC)
 
    % if bodies separating (i.e. if v_n>=0, then don't apply collision)
    sep_Idx = (v_n >= 0);               % (1 x NC)  (bool)
    sep_Idx = find(sep_Idx);            % indices   (int)
    v_up(:,sep_Idx) = v(:,sep_Idx);     % set to unchanged velocity v
end

function [Idx_co] = ClosestObjNode(x,x_co)
%{
function to find the closest object node 
to each particle/grid node inside the "collision body"
 
  NC = # colliding particles/grid nodes
  NObj = # object nodes

Input:  x       (NC x 3)    particle/grid node positions
        x_co    (NObj x 3)  object node positions  
  
Ouput:  Idx_co  (NC x 1)    indices of closest object node
%}
% Compute: 
    % displacement
    dx = x(:,1) - (x_co(:,1))';     % (NC x NObj)
    dy = x(:,2) - (x_co(:,2))';
    dz = x(:,3) - (x_co(:,3))';
    % distance
    D = ( dx.*dx + dy.*dy + dz.*dz ).^ 0.5;
    [~, Idx_co] = min(D,[],2);
end






                    