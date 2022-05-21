% mochi_pound_3d.m 
% Pounding Mochi 
% 3d Mochi model (coupled)

check_VE = 1;       % check when mochi is homogenous

d = 3;              % dimension
CFL = 0.25;
tf  = 10.0;
g = 1.0;            % gravitational constant
dtk = 0.0100;       % (dtk=0.0050, 200fps)(dtk=0.00625, 160fps)
                    % (dtk=0.0100, 100fps)(dtk=0.0167,  60fps) 
                    % (dtk=0.0333, 30fps) (dtk=0.0667,  15fps)

% material plot settings
h = [1.0 1.0 1.0];  % grid stepsize 
ppc  = 2;           % # particles per grid cell = ppc^d
psize = 25;         % for plotting particles
sqSz = 100;         % for plotting grid nodes

% Viscoelastic constants (Mochi)
md_p = 1.28;                                  % particle mass density
mu_min = 3;             mu_max = 6;           % shear modulus
La_min = 10*mu_min;     La_max = 10*mu_max;   % Lame's 1st parameter
mu_N_min = 1.0*10^-4;   mu_N_max = 1.0*10^-2; % viscous Newton constant
Wi_min = 0.75;          Wi_max = 5;           % Weissenberg number
Xi1 = 0.2; Xi2 = 0.5*10^-3; Xi3 = 0.2;        % "mochify parameters" 
Qup = @(Q,Jk) Qconst_up(Q,Jk,mu_min,mu_N_max,Wi_min,Xi1,Xi2,Xi3); % funct. handle

% ========================== MALLET SETTINGS =============================
% Mallet size (1.5,6)
r = 1.25;            % radius
bhalf = 4.0;         % half base   
maxMal = round((4*r*r)*(2*bhalf+1));       % max # mallet nodes

% Mallet Times and "Speed" controls
mochiWait = 0.0;    % wait for mochi to settle before pounding
rest = 0.25;        % time resting at top and bottom of swing
boost = 7.5;        % speed up swing (10 = 1 pound/sec)
wait = 0.75;        % wait to start swinging (mallet2)
                    % (boost,wait) = (3,2.5), (20,0.5), (15,0.75) (17,0.75)

% Rightward mallet (mallet1)
parabola1.path = 1; 
parabola1.a = -2/9; parabola1.h = 3; parabola1.k = 14;
parabola1.upper = [0, 8, 12]; parabola1.lower = [9, 8, 6];
parabola1.rest = rest; parabola1.boost = boost;
mallet1.r = r; 
mallet1.bhalf = bhalf;
mallet1.center = parabola1.upper'; 
mallet1.dir = [parabola1.path 0 0]';
mallet1.v_co = [0 0 0]';

% Leftward mallet (mallet2)
parabola2.path = -1; 
parabola2.a = -2/9; parabola2.h = 16; parabola2.k = 14;
parabola2.upper = [19, 11, 12]; parabola2.lower = [10, 11, 6];
parabola2.rest = rest; parabola2.boost = boost;
mallet2.r = r; 
mallet2.bhalf = bhalf;
mallet2.center = parabola2.upper'; 
mallet2.dir = [parabola2.path 0 0]';
mallet2.v_co = [0 0 0]';
                    
% ============================== MOCHI =================================== 
% Mochi initial values (sphere)
    c1x = 9.5;  c1y = 9.5; c1z = 5.5;  % center = (9.5,9.5,5.5)
    v1x = 0.0;  v1y = 0.0; v1z = 0.0;  % velocity = <0.0,0.0>
    r1in = 0.0; r1out = 4.0;           % inner, outer radius

% ======================= LAGRANGIAN PARTICLES ===========================
    % Mochi Particles 
    xc = (c1x-r1out) + h(1)/(2*ppc) : h(1)/ppc : (c1x+r1out) - h(1)/(2*ppc);  
    yc = (c1y-r1out) + h(2)/(2*ppc) : h(2)/ppc : (c1y+r1out) - h(2)/(2*ppc);
    zc = (c1z-r1out) + h(3)/(2*ppc) : h(3)/ppc : (c1z+r1out) - h(3)/(2*ppc);
    
    [Xi,Xj,Xk] = meshgrid(xc,yc,zc);  % c1: particle positions (in pink box)
        % Remove particles "outside" of circle c1 
    out = ( Xi(:)-c1x ).^2 + ( Xj(:)-c1y ).^2 + ( Xk(:)-c1z ).^2;
    out = [find( double(out > r1out^2) );       % outside particle indices
           find( double(out < r1in^2) )  ];
    Xi = Xi(:); Xj = Xj(:); Xk = Xk(:);         % reshape to col vectors
    Xi(out) = []; Xj(out) = []; Xk(out) = [];   % remove outside particles
        % Set initial particle positions
    Phi_k = [Xi, Xj, Xk];                       % current particle coords
    NP = size(Phi_k,1);
    
    % VE params struct for funct. handles 
    Qconst.mu = mu_max*ones(1,NP);          % decreases
    Qconst.La = La_max*ones(1,NP);          % decreases
    Qconst.mu_N = mu_N_min*ones(1,NP);      % increases
    Qconst.Wi = Wi_max*ones(1,NP);          % decreases
    
% ======================== EULERIAN MESH ================================= 
% h = [1.0, 1.0, 1.0];            % grid stepsize  
% h = [0.5, 0.5, 0.5];            % finer mesh

% grid struct
grid.h = h;
[grid,Xg,NG] = EulGrid(Phi_k,grid);
grid.height = Xg(:,3)-1;              % *** minus USU height = 1 ***

% ======================== MALLETS & USU =================================

% DEFAULT Mallet (centered at origin)
[Phi_mallet, n_mallet] = MalletLevelSet(r,bhalf,grid,3);    % ppc = 3

% Mallet Initial Values (t=0) 
[mallet1] = MalletUp(mallet1, parabola1, 0, dtk, Phi_mallet, n_mallet);
[mallet2] = MalletUp(mallet2, parabola2, 0, dtk, Phi_mallet, n_mallet);

% Usu (static)
usu.ctr = [c1x c1y 7]; usu.r = 6.0; usu.ht = 6.0;
[usu,cavIdx,topIdx] = UsuLevelSet(usu,grid,3);  % ppc = 3

% ========================= CHECK GRAPHS =================================
figure(1)
scatter3(Phi_k(:,1),Phi_k(:,2),Phi_k(:,3),psize,'filled','c'); hold on; % mochi
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','b'); hold on; % grid
M1C = [mallet1.center]; M2C = [mallet2.center]; % mallet centers
scatter3([M1C(1,end); M2C(1,end)],[M1C(2,end); M2C(2,end)],[M1C(3,end); M2C(3,end)],'k','filled'); hold on;
scatter3(mallet1.x_co(:,1),mallet1.x_co(:,2),mallet1.x_co(:,3),'ro'); hold on; % mallet1
scatter3(mallet2.x_co(:,1),mallet2.x_co(:,2),mallet2.x_co(:,3),'bo'); hold on; % mallet2
scatter3(usu.x_co(:,1),usu.x_co(:,2),usu.x_co(:,3),1,'ko'); hold off; % usu
xlabel('x'); ylabel('y'); zlabel('z');
xlim([-4 23]); ylim([0 19]); zlim([0 17]);

% ============================ Step k = 0 ================================
k = 0;
t = 0;      % time at step k

% Compute Particle Volumes (Tot_Volume/NP)
    Vol = prod(h)*(1/ppc^d)*ones(NP,1); 
    
% Set Particle Mass 
    mass_p = md_p*ones(NP,1).*Vol;  
    
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
    % funct. handle: grid, Vol, mass_p, & g are constant over time
    RealEngy = @(J,B_El,Uk,w,Xp,Qconst) RealEngy_comp(J,B_El,Uk,w,Xp,grid,Vol,mass_p,Qconst,g);
    u = P2G_vec(grid,Phi_k,Uk,[0 0 0]); w = u; % so we can compute init. Engy
    
    % ANIMATION
    T = t; 
    PHX = Phi_k(:,1); 
    PHY = Phi_k(:,2);
    PHZ = Phi_k(:,3);
    JPK = det_3d(Fk)';
    M1X = mallet1.x_co(:,1);    M2X = mallet2.x_co(:,1);
    M1Y = mallet1.x_co(:,2);    M2Y = mallet2.x_co(:,2);
    M1Z = mallet1.x_co(:,3);    M2Z = mallet2.x_co(:,3);
    M1C = mallet1.center;       M2C = mallet2.center;

% ============================= Step k =================================== 
tic
while t < tf
    % "Actual" Energy over time
    EvT(k+1,:) = [t,RealEngy(Jk,Bk,Uk,w,Phi_k,Qconst)]; % [t,E,Q,K,Q_El,Q_N,G]
    
    %  APIC CFL (non-uniform time steps)
%     dtk = CFL*h(1)/( max(sum(abs(Uk),1) + (6/h(1))*sqrt(d)*sqrt(sum(Bm.^2,1))) );
%     dtk = min(dtk,tf-t);
    
    % MALLETS
    if t >= mochiWait           % Mallet 1      
        [mallet1] = MalletUp(mallet1,parabola1,t-mochiWait,dtk,Phi_mallet,n_mallet);
    end
    if t >= wait + mochiWait    % Mallet 2 (wait a bit after Mallet 1)               
        [mallet2] = MalletUp(mallet2,parabola2,t-(mochiWait+wait),dtk,Phi_mallet,n_mallet); 
    end
    
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
    [Idx1] = CollisionDetect_Mallet(Xg,mallet1,grid);
    [Idx2] = CollisionDetect_Mallet(Xg,mallet2,grid);
    [Idx3] = CollisionDetect_Usu(Xg,usu,grid);
    [v_up1] = CollisionHandling(Idx1,Xg(Idx1,:),w(:,Idx1),mallet1);
    [v_up2] = CollisionHandling(Idx2,Xg(Idx2,:),w(:,Idx2),mallet2);
    [v_up3] = CollisionHandling(Idx3,Xg(Idx3,:),w(:,Idx3),usu);
    w(:,Idx1) = v_up1; w(:,Idx2) = v_up2; w(:,Idx3) = v_up3;
    
    % APIC TRANSFER grid to particle (advection) 
    xup = Xg' + dtk*w;
    [Phi_up,Vk,Bm] = G2P_APIC(grid,Phi_k,xup,w); 
    
    % COLLISION HANDLING: particle velocities
    [Idx1] = CollisionDetect_Mallet(Phi_up,mallet1,grid);
    [Idx2] = CollisionDetect_Mallet(Phi_up,mallet2,grid);
    [Idx3] = CollisionDetect_Usu(Phi_up,usu,grid);
    [v_up1] = CollisionHandling(Idx1,Phi_up(Idx1,:),Vk(:,Idx1),mallet1);
    [v_up2] = CollisionHandling(Idx2,Phi_up(Idx2,:),Vk(:,Idx2),mallet2);
    [v_up3] = CollisionHandling(Idx3,Phi_up(Idx3,:),Vk(:,Idx3),usu);
    Vk(:,Idx1) = v_up1; Vk(:,Idx2) = v_up2; Vk(:,Idx3) = v_up3;
    
    % Particle Advection (collision corrected Lag velocity)  
    Phi_up = Phi_k + dtk*Vk';

    % GRAPH
    scatter3(Phi_up(:,1),Phi_up(:,2),Phi_up(:,3),psize,'filled','co'); hold on; % mochi
    scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on; % grid
    scatter3([M1C(1,end); M2C(1,end)],[M1C(2,end); M2C(2,end)],[M1C(3,end); M2C(3,end)],'k','filled'); hold on; % mallet centers
    scatter3(mallet1.x_co(:,1),mallet1.x_co(:,2),mallet1.x_co(:,3),1,'ro'); hold on; % mallet1
    scatter3(mallet2.x_co(:,1),mallet2.x_co(:,2),mallet2.x_co(:,3),1,'bo'); hold on; % mallet2
    scatter3(usu.x_co(:,1),usu.x_co(:,2),usu.x_co(:,3),1,'ko'); hold off; % usu
    title(num2str(kb),num2str(t)); xlabel(num2str(k)); ylabel(num2str(dtk));
    xlim([-4 23]); ylim([0 19]); zlim([0 17]);
    hold off; pause(1e-10);
    
    % ANIMATION
    T = [T,t]; 
    PHX = [PHX, Phi_up(:,1)]; 
    PHY = [PHY, Phi_up(:,2)]; 
    PHZ = [PHZ, Phi_up(:,3)];
    JPK = [JPK, det_3d(Fk)'];
    M1X = [M1X, mallet1.x_co(:,1)]; M2X = [M2X, mallet2.x_co(:,1)];
    M1Y = [M1Y, mallet1.x_co(:,2)]; M2Y = [M2Y, mallet2.x_co(:,2)];
    M1Z = [M1Z, mallet1.x_co(:,3)]; M2Z = [M2Z, mallet2.x_co(:,3)];
    M1C = [M1C, mallet1.center];    M2C = [M2C, mallet2.center];
    
    % Update vars for next time step 
    k = k + 1;
    t = t + dtk;
    Fk = F_comp(grid,Phi_k,dtk,w,Fk);   % 9xNP
    Jk = det_3d(Fk);                    % 1xNP
    Bob_k = Bob_comp(Bob_k,dtk,Qconst,grid,Phi_k,w);    % 9xNP
    Bk = B_comp(Fk,Bob_k);              % 9xNP
    Phi_k = Phi_up;                     % NPx3
    Uk = Vk;                            % 3xNP
    Qconst = Qup(Qconst,Jk);            % VE params
    
    % Update Eulerian grid
    [grid,Xg,NG] = EulGrid(Phi_k,grid);

    % Output: t & k (if Mochi is homogenous)
    [check_VE] = MochiHomog(check_VE,Qconst,mu_min,mu_N_max,Wi_min,t,k);
end
time = toc;

% "Actual" Energy over time 
EvT(k+1,:) = [t,RealEngy(Jk,Bk,Uk,w,Phi_k,Qconst)]; % [t,E,Q,K,Q_El,Q_N,G]

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

% Write J=det(Fp) particle data to text file
writematrix(JPK,'Jpk.dat','Delimiter',';');  % (NP x k+1)

% Write Usu Node Positions to text file
writematrix(usu.x_co,'USU.dat','Delimiter',';');

% Write Mallet Node Positions to text file
writematrix(M1X,'M1X.dat','Delimiter',';'); 
writematrix(M1Y,'M1Y.dat','Delimiter',';');
writematrix(M1Z,'M1Z.dat','Delimiter',';');
writematrix(M2X,'M2X.dat','Delimiter',';');
writematrix(M2Y,'M2Y.dat','Delimiter',';');
writematrix(M2Z,'M2Z.dat','Delimiter',';');

% Write Mallet centers to text file
writematrix(M1C,'M1C.dat','Delimiter',';');
writematrix(M2C,'M2C.dat','Delimiter',';'); 

% -------------------------- ANIMATED SIM --------------------------------
% SAVE Workspace and IMPORT into mochi_pound_3d_ani.m  
% in order to render simulation as gif/mp4

% ++++++++++++++++++++++++++++ FUNCTIONS +++++++++++++++++++++++++++++++++
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

function [Q] = Qconst_up(Q,Jk,mu_min,mu_N_max,Wi_min,Xi1,Xi2,Xi3)
%{
Q struct of VE params:
    mu                  1xNP
    La                  1xNP
    mu_N                1xNP
    Wi                  1xNP
Jk                      1xNP
%}
NP = length(Jk);

Q.mu = max( mu_min*ones(1,NP), Q.mu - Xi1*(abs(Jk-1)) );
Q.La = Q.mu*10;
Q.mu_N = min( mu_N_max*ones(1,NP), Q.mu_N + Xi2*(abs(Jk-1)) );
Q.Wi = max( Wi_min*ones(1,NP), Q.Wi - Xi3*(abs(Jk-1)) );

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
    Q_N = ( mu_N.* dot(S,S,1) )*Vol_k;% Newtonian viscous potential energy
    G = g*dot(mass_p, Xp(:,3)-1);     % Gravitational potential energy *** minus USU height = 1 ***
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

% ========================== MOCHI FUNCTIONS ==============================
function [Phi] = MochiFlip(Phi,h,t)
% Rotate particles 180deg in xz-plane (3d) or xy-plane (2d) about CoM 

% if t = 0, don't flip
if (t==0)
   return; 
end

% Find center of mass in xz-plane (3d) or xy-plane (2d)
CoM = mean( [Phi(:,1), Phi(:,end)], 1 );

% Compute vertical "lift" (so mochi is within usu cavity)
lift = h(end); %max(Phi(:,end)) - CoM(2); 

% Set new x & z-coords (3d) or x & y-coords (2d)
Phi(:,1) = -Phi(:,1) + 2*CoM(1);
Phi(:,end) = -Phi(:,end) + 2*CoM(2) + lift;  
end

function [Phi_up] = MochiRotate(Phi,usuCtr,t)
% Rotate particles 90deg in xy-plane (3d) about usuCtr 
% z-coord stays the same

% if t = 0, don't rotate
if (t==0)
   return; 
end

% Create temp array to hold new Particle Positions
Phi_up = Phi;

% Set new x & y-coords (3d)
Phi_up(:,1) = -Phi(:,2) + ( usuCtr(2) + usuCtr(1) );
Phi_up(:,2) =  Phi(:,1) + ( usuCtr(2) - usuCtr(1) );
end

function [check_VE] = MochiHomog(check_VE,Q,mu_min,mu_N_max,Wi_min,t,k)
% function to check if mochi particles are homogenous 

% if Mochi is ALREADY homogenous (i.e. check_VE = 0)
    if( check_VE == 0 )  % don't check VE params 
       return;           % exit function
    end

% else: check if homogeneous
    mu = Q.mu;      ep1 = 0.1;
    mu_N = Q.mu_N;  ep2 = 0.1*10^-3;
    Wi = Q.Wi;      ep3 = 0.1;

    % if homogenous
    if( max(mu) < mu_min + ep1 && min(mu_N) > mu_N_max - ep2 && max(Wi) < Wi_min + ep3 )
        disp('t= '); disp(t);   % print time
        disp('k= '); disp(k);   % print iter
        check_VE = 0;           % don't check VE params anymore
    end
end

% ================= MALLET FUNCTIONS (KINEMATIC MALLET) ==================
% Mallet Centers from parabola path
function [centerUp] = malletCenter(parabola,t)
%{
function that returns the center coordinate of the mallet at time t,
where mallet's path is given by the parabola struct

*** keeping  y-coord static for all t ***
parabola (struct):
    f(x) = a(x-h)^2 + k
    upper = (x,y,height)
    lower = (x,y,height)
    path: 1 or -1 (1: rightward, -1: leftward) (during descent)
    rest: time duration at rest (top and bottom)
    boost = speed boost

Output: centerUP    (3 x 1) 
%}

% Unpack parabola struct
a = parabola.a; h = parabola.h; k = parabola.k;
upper = parabola.upper; lower = parabola.lower;
path = parabola.path; rest = parabola.rest; boost = parabola.boost;

% y-coord static
y = upper(2);

xdom = lower(1) - upper(1); % x domain

tau_max = 2*(abs(xdom) + rest);
tau = boost*mod(t,tau_max);         
tau = mod(tau,tau_max);     % keep positive

% rightward path
if path > 0
    % 1st pause
    if tau <= rest
        xUp = upper(1);
        f = upper(end);
    % descent
    elseif tau <= rest + xdom
        xUp = upper(1) + (tau - rest);
        f = a*(xUp - h)^2 + k;
    % 2nd pause    
    elseif tau <= 2*rest + xdom
        xUp = lower(1);
        f = lower(end);
    % ascent
    else
        xUp = lower(1) - (tau - 2*rest - xdom);
        f = a*(xUp - h)^2 + k; 
    end
end

% leftward path
if path < 0
    % 1st pause
    if tau <= rest
        xUp = upper(1);
        f = upper(end);
    % descent
    elseif tau <= rest + abs(xdom)
        xUp = upper(1) - (tau - rest);
        f = a*(xUp - h)^2 + k;
    % 2nd pause
    elseif tau <= 2*rest + abs(xdom)
        xUp = lower(1);
        f = lower(end);
    % ascent
    else
        xUp = lower(1) + (tau - 2*rest - abs(xdom));
        f = a*(xUp - h)^2 + k; 
    end    
end

centerUp = [xUp,y,f]';
end

% Initial Level Set
function [Phi_mallet,n,frIdx,bkIdx,sdIdx] = MalletLevelSet(r,bhalf,grid,ppc)
%{ 
function to compute Level Set nodes and normals for the mallet 
center:(0,0,0), front & back: yz-plane, horizontal in xz-plane

    NObj = # object nodes

Input:  r       mallet radius
        bhalf   half of mallet lenght
        grid    Eulerian grid struct   
        ppc     particles per cell (for mallet)

Output: Phi_mallet          (NObj x 3) node positions
        n                   (3 x NObj) normal vectors
        frIdx,bkIdx,sdIdx   front, back, & side node indices (debugging)
%}

% unpack grid struct
h = grid.h;

% Rectangle 
yc = -round(r+0.5)+ h(2)/(2*ppc): h(2)/ppc : round(r+0.5) - h(2)/(2*ppc);
zc = -round(r+0.5)+ h(3)/(2*ppc): h(3)/ppc : round(r+0.5) - h(3)/(2*ppc);
[Xj,Xk] = meshgrid(yc,zc);

% Circle (center: origin, radius: r)
    % Remove particles "outside" of circle (center: origin, radius r) 
    out = Xj(:).^2 + Xk(:).^2;
    out = find( double(out > (r+0.05)^2) );     % outside particle indices 
    Xj_cir = Xj(:); Xk_cir = Xk(:);             % reshape to col vectors
    Xj_cir(out) = []; Xk_cir(out) = [];         % remove outside particles
    N_cir = length(Xj_cir);                     % # nodes in circle
    
% Annulus 
    Xj_ann = Xj_cir; Xk_ann = Xk_cir;
    out = Xj_ann(:).^2 + Xk_ann(:).^2;
    out = find( double(out < (r-0.20)^2) );     % outside particle indices 
    Xj_ann(out) = []; Xk_ann(out) = [];         % remove outside particles
    N_ann = length(Xj_ann);                     % # nodes in annulus
    n_ann = [zeros(1,N_ann);                    % annulus (side) normals
             Xj_ann';                           %       (3 x N_ann)
             Xk_ann'];  
    n_ann_norm = vecnorm(n_ann,2,1);
    n_ann = bsxfun(@rdivide,n_ann,n_ann_norm);  % unit length normals

% Mallet Front & Back (disks)
    % Nodes
    xfr = -bhalf + h(1)/(2*ppc);
    xbk =  bhalf - h(1)/(2*ppc);
    Phi_mallet = [xfr*ones(N_cir,1),Xj_cir, Xk_cir;     % front
                  xbk*ones(N_cir,1),Xj_cir, Xk_cir];    % back
    % Indices
    frIdx = 1: 1: N_cir; 
    bkIdx = N_cir + 1: 1: 2*N_cir;
    % Normals
    %    front                      back
    n = [repmat([-1 0 0]',1,N_cir), repmat([1 0 0]',1,N_cir)];

% Mallet Side (extruded annulus)
    % Nodes
    xsd = -bhalf + 1.5*(h(1)/ppc): h(1)/ppc: bhalf - 1.5*(h(1)/ppc);
    N_extr = length(xsd);       % # annuli to extrude
    xsd = repmat(xsd,N_ann,1);
    xsd = xsd(:);
    Phi_mallet = [Phi_mallet;   % front & back
             xsd, repmat(Xj_ann,N_extr,1), repmat(Xk_ann,N_extr,1)]; % side
    % Indices
    sdIdx = 2*N_cir + 1: 1: size(Phi_mallet,1);
    
    % Normals
    %    fr & bk    side
    n = [n,         repmat(n_ann,1,N_extr)];
    
    figure(4);
    scatter3(Phi_mallet(frIdx,1),Phi_mallet(frIdx,2),Phi_mallet(frIdx,3),'ro'); hold on;
    scatter3(Phi_mallet(bkIdx,1),Phi_mallet(bkIdx,2),Phi_mallet(bkIdx,3),'md'); hold on;
    scatter3(Phi_mallet(sdIdx,1),Phi_mallet(sdIdx,2),Phi_mallet(sdIdx,3),'gs'); hold off;
    title('Mallet Level Set');
end

% Update & Transform Mallet 
function [mallet] = MalletUp(mallet,parabola,t,dt,Phi_mallet,n_mallet)
%{
function to update mallet struct 
(update center, dir, translate nodes, compute velocity, rotate normals)

mallet: center  (3 x 1) 
        dir     (3 x 1) *** must be unit length ***
        x_co    (NObj x 3)
        v_co    (3 x 1)
        n_co    (3 x NObj)

Input:  mallet (struct)
        parabola (struct)
        t
        dt
        Phi_mallet          (NObj x 3)  DEFAULT mallet centered at (0,0,0)
        n_mallet            (3 x NObj)  DEFAULT mallet normals

Output: mallet (struct)
%}
% Unpack mallet struct
    center = mallet.center;
    dir = mallet.dir;    


% Update mallet center
    oldCent = center;
    centerUp = malletCenter(parabola,t);   % (3 x 1)  

% Compute: direction and velocity
    v_co = (centerUp - oldCent)/dt;
    % if mallet stationary (i.e. v_co = 0), don't update dir    
    if ( norm(v_co,2) > eps )
        dir = v_co/( norm(v_co,2) );     
    end
    
% Rotate mallet & normals about origin (in xz-plane only) 
    % *** dir must be unit length ***
    Phi_mallet = RotateXZ(Phi_mallet',dir)'; 
    n_co = RotateXZ(n_mallet,dir);        

% Translate mallet
    Phi_mallet = bsxfun(@plus,Phi_mallet,centerUp'); 

% Update struct 
    mallet.center = centerUp;
    mallet.dir = dir;
    mallet.x_co = Phi_mallet;
    mallet.v_co = v_co;
    mallet.n_co = n_co;
end

% 2D Rotation function 
function [vectRot,rotAng] = RotateXZ(vect,dir)
%{
function to rotate array of 3d vectors in CCW xz-plane, about origin
"Rotation Angle" is with respect to e1 = <1 0 0>
where dir is ONLY in xz-plane 
    s.t. cross(dir,e1)= c*e2, constant c in [-1,1]
    i.e. vector parallel to y-axis

    theta   = acos( dot(dir,e1) ), in [0,pi]
    phi     = asin( c )          , in [-0.5pi, 0.5pi]

Input:  vect    (3 x NV)
        dir     (3 x 1)     *** assuming unit length ***

Output: vectRot (3 x NV)
        rotAng  rotation angle in radians 
%}

% Compute: theta and c
    e1 = [1 0 0]';
    theta = acos(dot(dir,e1));
    c = cross(dir,e1);
    c = c(2);

% Find "Rotation Angle" between dir and e1 (in xz-plne)
    % if dir is in   Q3 or Q4 (i.e. c<0),  and  rot. angle is 2pi-theta 
    if (c < 0)
        rotAng = 2*pi - theta;
        
    % else dir is in Q1 or Q2 (i.e. c>=0), then rot. angle is theta
    else
        rotAng = theta;
    end

% Apply rotation
    cR = cos(rotAng); 
    sR = sin(rotAng);

    vectRot = vect;     % 3 x NV
    vectRot(1,:) = cR*vect(1,:) - sR*vect(3,:);
    % y-coord stays the same 
    vectRot(3,:) = sR*vect(1,:) + cR*vect(3,:);
end

% Mallet Collision Detection 
function [Idx] = CollisionDetect_Mallet(x,mallet,grid)
%{
Collision Detection function of Mallet (cross section = circle)
*** Movement only in the xz-plane, so y-coords are const ***

Input:  x       (NPoG x 3) particle/grid node positions
        mallet  (struct)
        grid    (struct)

mallet: center  (3 x 1) 
        dir     (3 x 1) *** should be unit length ***
        x_co    (NObj x 3)
        v_co    (3 x 1)
        n_co    (3 x NObj)
        
Output: Idx (NC x 1)    indices of particle/grid nodes "inside" mallet
%}
% Unpack mallet struct
    center = mallet.center;
    dir = mallet.dir;
    r = mallet.r;               % *** scaled with grid spacing ***
    bhalf = mallet.bhalf;       % *** scaled with grid spacing ***
  
% Format (rescale) variables 
    r     = r*max(grid.h);        
    bhalf = bhalf*max(grid.h);     
    
% Compute: domain of mallet in y-axis    
    y_min = min(mallet.x_co(:,2));
    y_max = max(mallet.x_co(:,2));
    
% Compute: cyllinder values
    base = 2*bhalf;
    tail = center - bhalf*dir;      % (3 x 1) blue point
    tip  = center + bhalf*dir;      % (3 x 1) green point

% for each particle/grid node 
    a = bsxfun(@minus,x',tail);         % (3 x NPoG) blue vector
    b = bsxfun(@minus,x',tip);          % (3 x NPoG) green vector
    abcross = cross(a,b);
    height = vecnorm(abcross,2,1)/base; % (1 x NPoG)

    c = bsxfun(@minus,x',center);               % (3 x NPoG)
    proj = abs(sum(bsxfun(@times,c,dir),1));    % (1 x NPoG) |dot(c,dir)|
    
% Check if inside mallet
    valid = (height' <= r).*(proj' <= bhalf);               % (NPoG x 1)
    valid = valid.*(x(:,2) >= y_min).*(x(:,2) <= y_max);    % (NPoG x 1)
    Idx = find(valid);                                      % (NC x 1)
end

% ==================== USU FUNCTIONS (STATIC) ======================
% Usu Level Set 
function [usu,cavIdx,topIdx] = UsuLevelSet(usu,grid,ppc)
%{ 
function to compute Level Set nodes and normals for the usu
(no nodes for the usu sides)

    NObj = # object nodes

Input:  usu.ctr,r,ht (incomplete struct)
        grid    Eulerian grid struct   
        ppc     particles per cell (for usu)

Output: usu (complete struct)
        cavIdx,topIdx   cavity & top node indices (debugging)

usu:    ctr         (1 x 3) center of usu cavity
        r           radius of usu cavity
        ht          height of usu (i.e. z-coord of usu top)
        x_co        (NObj x 3)  node positions
        v_co        (3 x 1)     node velocities
        n_co        (3 x NObj)  normal vectors
%}

% unpack grid struct
    h = grid.h;
    
% unpack usu struct
    usuCtr = usu.ctr;
    usuRad = usu.r;
    usuHt  = usu.ht;

% Cavity 
    xc = round(usuCtr(1)-usuRad)-1 + h(1)/(2*ppc) : h(1)/ppc : round(usuCtr(1)+usuRad)+1;  
    yc = round(usuCtr(2)-usuRad)-1 + h(2)/(2*ppc) : h(2)/ppc : round(usuCtr(2)+usuRad)+1;  
    zc = round(usuCtr(3)-usuRad)-1 + h(3)/(2*ppc) : h(3)/ppc : round(usuCtr(3)+usuRad)+1;  
    [Xi,Xj,Xk] = meshgrid(xc,yc,zc);  % box
    out = ( Xi(:)-usuCtr(1) ).^2 + ( Xj(:)-usuCtr(2) ).^2 + ( Xk(:)-usuCtr(3) ).^2;   % sq dist from center
    out = [find( double(out > (usuRad+0.05)^2) );   % outside shell indices 
           find( double(out < (usuRad-0.15)^2) )];  
    Xi = Xi(:); Xj = Xj(:); Xk = Xk(:);         % reshape to col vectors
    Xi(out) = []; Xj(out) = []; Xk(out) = [];   % remove out nodes
    above = find( double(Xk > usuHt) );         % above usuHt
    Xi(above)=[]; Xj(above)=[]; Xk(above)=[];   % remove above nodes

    % Nodes
    Phi_usu = [Xi, Xj, Xk];     % cavity nodes
    N_cav = size(Phi_usu,1);    
    % Indices
    cavIdx = 1: 1: N_cav;
    % Normals 
    n_cav = bsxfun(@minus,Phi_usu,usuCtr)';        % shift center to origin
    n_cav_norm = vecnorm(n_cav,2,1);
    n_cav = bsxfun(@rdivide, n_cav, n_cav_norm);   % unit length normals
    n_cav = - n_cav;                               % "inward" 

% Top
    a1 = grid.O(1)-2; a2 = max(grid.x)+2; 
    b1 = grid.O(2)-2; b2 = max(grid.y)+2;
    xc = a1 + h(1)/(2*ppc) : h(1)/ppc : a2 - h(1)/(2*ppc);  
    yc = b1 + h(2)/(2*ppc) : h(2)/ppc : b2 - h(2)/(2*ppc);    
    [Xi,Xj] = meshgrid(xc,yc);  % rectangle
    hRad = usuRad + 0.1*h(1);   
    hole = ( Xi(:)-usuCtr(1) ).^2 + ( Xj(:)-usuCtr(2) ).^2; 
    hole = find( double(hole < hRad^2) );
    Xi = Xi(:); Xj = Xj(:);         % reshape to col vectors
    Xi(hole) = []; Xj(hole) = [];   % remove hole nodes
    
    % Nodes
    N_top = length(Xi);                         % # top nodes
    Phi_usu = [Phi_usu;                         % cavity  nodes
               Xi, Xj, usuHt*ones(N_top,1)];    % top nodes
    % Indices
    topIdx = N_cav+1: 1: N_cav+N_top;
    % Normals
    %    cav    top
    n = [n_cav, repmat([0 0 1]',1,N_top)];
    
% "Update" usu struct (to include new vars)
    usu.x_co = Phi_usu;
    usu.v_co = [0 0 0]';
    usu.n_co = n;
end

% Usu Collision Detection  
function [Idx] = CollisionDetect_Usu(x,usu,grid)
%{
Collision Detection function of Usu (static)

Input:  x       (NPoG x 3) particle/grid node positions
        usu     (struct)
        grid    (struct)

usu:    ctr         (1 x 3) center of usu cavity
        r           radius of usu cavity
        ht          height of usu (i.e. z-coord of usu top)
        x_co        (NObj x 3)  node positions
        v_co        (3 x 1)     node velocities
        n_co        (3 x NObj)  normal vectors
        
Output: Idx (NC x 1)    indices of particle/grid nodes "inside" usu
%}
        
% Unpack usu struct
    ctr = usu.ctr;
    r   = usu.r;        
    ht  = usu.ht;
    xmin = min(usu.x_co(:,1));
    xmax = max(usu.x_co(:,1));
    ymin = min(usu.x_co(:,2));
    ymax = max(usu.x_co(:,2));
    
% Compute: distance from usu center (for each particle/grid node)
    duc = bsxfun(@minus,x,ctr);         % (NPoG x 3)
    duc = sqrt( sum( duc.*duc, 2 ) );   % (NPoG x 1)

% Check if inside usu    
    in = (duc > r) .* (x(:,3) < ht);
    in = in .* (duc < r + 1.5*max(grid.h));  % since too many grid nodes
    in = in .* (x(:,1) > xmin) .* (x(:,1) < xmax);
    in = in .* (x(:,2) > ymin) .* (x(:,2) < ymax);
    Idx = find(in);
end

% ======================= COLLISION FUNCTIONS ============================
function [v_up] = CollisionHandling(Idx,x,v,object)
%{
% function to update the velocities of particles or grid nodes 
  "colliding" with the object (that undergoes "rigid" transformations)

  NC = # colliding particles/grid nodes
  NObj = # object nodes

Input:  x       (NC x 3)    particle/grid node positions
        v       (3 x NC)    particle/grid node velocities
        object  (struct)

object: v_co    (3 x 1)     object velocity (same for all nodes)
        x_co    (NObj x 3)  object node positions
        n_co    (3 x NObj)  object node normals (*ALL nodes are bdry nodes)

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
    
% Unpack object struct
    v_co = object.v_co;
    x_co = object.x_co;
    n_co = object.n_co;
    
% friction coeff
    mu = 0.3;  % sliding friction coeff. (wood)
%     mu = 0.5;   % static friction coeff. (wood)

% Compute: 
    v_rel = bsxfun(@minus, v, v_co);    % relative velocity (3 x NC)
    
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

          


