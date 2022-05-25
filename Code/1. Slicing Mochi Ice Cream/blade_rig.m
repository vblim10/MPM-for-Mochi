% blade_rig.m
% Script to direct the positioning & movement of the kinematic blades
% before running simulation with mochi ice cream (mochi_midairSlice.m)

% *** IMPORT WORKSPACE *** wkspc_mochi_ctr_xy_8_8.mat ***

d = 3;              % dimension
tf  = 2.3;
dtk = 0.0166;       % (dtk=0.0166, 60 fps)
                    % (dtk=0.0333, 30 fps) 
                    % (dtk=0.0667, 15 fps)

% Plot settings
psize = 25;         % for plotting particles
sqSz = 100;         % for plotting grid nodes

% ======================= LAGRANGIAN PARTICLES ===========================
% *** IMPORT WORKSPACE *** (Phi_up, Idx_cream, Idx_mochi)
% Particle postions for Ice cream & Mochi
Phi_k = Phi_up + [0 0 5];   % Shift upward
NP = size(Phi_k,1);
NP_cream = length(Idx_cream);
NP_mochi = length(Idx_mochi);

% ======================== EULERIAN MESH ================================= 
h = [1.0, 1.0, 1.0];            % grid stepsize  
% h = [0.5, 0.5, 0.5];            % finer mesh

% grid struct
grid.h = h;
[grid,Xg,NG] = EulGrid(Phi_k,grid);
    
% ======================== KINEMATIC BLADE =============================== 
% Blade constants
hb = 0.32;  % blade node spacing & thickness
a = 7;      % ellipse long "radius" (scaled with hb)

% Default Blade
[blade1,Phi_blade,n_co] = BladeLevelSet(a,hb,grid);
    
% Blade 1
blade1.d = [-1 -1]./sqrt(2);
blade1.o = [15.5 16.5];
blade1.x_co = RotateXZ(Phi_blade', blade1.d)'+ [blade1.o(1) 0 blade1.o(2)];    
blade1.v_co = [-11; 0; -11];    % init. vel = <-11,0,-11>   (3 x 1)
blade1.n_co = RotateXZ(n_co, blade1.d); 
blade1.a_co = [0; 0; 0];        % (3 x 1)

% Blade 2
blade2 = blade1;
blade2.d = [1 -1]./sqrt(2);
blade2.o = [-3 20];
blade2.x_co = RotateXZ(Phi_blade', blade2.d)' + [blade2.o(1) 0 blade2.o(2)];    
blade2.v_co = [11; 0; -11];      % init. vel = <11,0,-11>    (3 x 1)
blade2.n_co = RotateXZ(n_co, blade2.d); 
blade2.a_co = [0; 0; 0];        % (3 x 1)  

% Blade 3
blade3 = blade1;
blade3.d = [-1 0];
blade3.o = [32 8];
blade3.x_co = RotateXZ(Phi_blade', blade3.d)' + [blade3.o(1) 0 blade3.o(2)];    
blade3.v_co = [-15.5; 0; 0];    % init. vel = <-15.5,0,0>   (3 x 1)
blade3.n_co = RotateXZ(n_co, blade3.d); 
blade3.a_co = [0; 0; 0];        % (3 x 1)

% =========================== CHECK GRAPH ================================
figure(1)
scatter3(Phi_k(Idx_cream,1),Phi_k(Idx_cream,2),Phi_k(Idx_cream,3),psize,'filled','m'); hold on; % ice cream
scatter3(Phi_k(Idx_mochi,1),Phi_k(Idx_mochi,2),Phi_k(Idx_mochi,3),psize,'filled','g'); hold on; % mochi
scatter3(blade1.x_co(:,1),blade1.x_co(:,2),blade1.x_co(:,3),psize,'filled','k'); hold on; % blade1
scatter3(blade2.x_co(:,1),blade2.x_co(:,2),blade2.x_co(:,3),psize,'filled','k'); hold on; % blade2
scatter3(blade3.x_co(:,1),blade3.x_co(:,2),blade3.x_co(:,3),psize,'filled','k'); hold on; % blade3
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on;
az = 0.0447; el = 0.0677;   % side view (xz-plane)
view([az,el]);
hold off;

% Blade Rigging
t = 0.0; k = 0;
while t < tf

    % Plot
    figure(1);
    scatter3(Phi_k(Idx_cream,1),Phi_k(Idx_cream,2),Phi_k(Idx_cream,3),psize,'filled','m'); hold on; % ice cream
    scatter3(Phi_k(Idx_mochi,1),Phi_k(Idx_mochi,2),Phi_k(Idx_mochi,3),psize,'filled','g'); hold on; % mochi
    scatter3(blade1.x_co(:,1),blade1.x_co(:,2),blade1.x_co(:,3),psize,'filled','k'); hold on; % blade1 
    scatter3(blade2.x_co(:,1),blade2.x_co(:,2),blade2.x_co(:,3),psize,'filled','k'); hold on; % blade2
    scatter3(blade3.x_co(:,1),blade3.x_co(:,2),blade3.x_co(:,3),psize,'filled','k'); hold on; % blade3
    scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on;
    title(['$t_{k}=$',' ',num2str(t)],'Interpreter','latex'); 
    xlabel(num2str(k)); ylabel(num2str(dtk));
    xlim([-5 30]); ylim([0 16]); zlim([2 20]);
    az = 0.0447; el = 0.0677;   % side view (xz-plane)
    view([az,el]);  
    hold off; pause(1e-10);
    
    % update values
    [blade1] = BladeUp(blade1,dtk);
    [blade2] = BladeUp(blade2,dtk);
    [blade3] = BladeUp(blade3,dtk); 
    k = k+1;
    t = t + dtk;
    
end

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

% =================== BLADE FUNCTIONS (KINEMATIC BLADE) ==================
% Initial Level Set
function [blade,Phi_blade,n_co] = BladeLevelSet(a,hb,grid)
%{
function to compute Level Set nodes and normals for the blade 
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
% y-coords of nodes (constant)
blade.y = [min(grid.y)+0.5*hb: hb: max(grid.y)-0.5*hb]';    
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
blade.o = blade.o + dt*[v_co(1) v_co(3)]; % x,z-coord of ellipse center (1 x 2)

% Update node velocities
v_co = v_co + dt*a_co;

blade.v_co = v_co;
blade.x_co = x_co;
blade.n_co = n_co;   % stay same (no rotation)
end

function [Idx] = CollisionDetect_Blade(x,blade)
%{ 
Collision Detection function of lenticular blade (cross section=ellipse)
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

Input:  vect        (3 x NV)
        dir         (3 x 1)     *** assuming unit length ***

Output: vectRot     (3 x NV)
        rotAng      rotation angle in radians 
%}
% include y-dim for direction
dir = [dir(1) 0 dir(2)];

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
    
% Add pi/2 for blade *** (default blade is along [x z]=[0 -1]) ***
rotAng = rotAng + 0.5*pi;

% Apply rotation
    cR = cos(rotAng); 
    sR = sin(rotAng);

    vectRot = vect;     % 3 x NV
    vectRot(1,:) = cR*vect(1,:) - sR*vect(3,:);
    % y-coord stays the same 
    vectRot(3,:) = sR*vect(1,:) + cR*vect(3,:);
      
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





