% mallet_rig_3d.m
% Script to direct the positioning & movement 
% of the kinematic mallets & static usu
% before running simulation with mochi

tf = 10; 
time = linspace(0,tf,126); 
dt = time(2)-time(1);
k = 1;
t = time(k);

mochiWait = 0.0;    % wait for mochi to settle before pounding
rest = 0.25;        % time resting at top and bottom of swing
boost = 7.5;        % speed up swing (10 = 1 pound/sec)
wait = 0.75;        % wait to start swinging (mallet2)
                    % (boost,wait) = (3,2.5), (20,0.5), (15,0.75) (17,0.75)
         
% mallet size (r=1.5, bhalf=4.5, h=0.5)
r = 1.25;           % radius 
bhalf = 4.0;        % half base  

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

% MPM Eulerian grid 
h = [1.0, 1.0, 1.0];            % grid stepsize    
% h = [0.5, 0.5, 0.5];            % finer mesh
x = 2.0:h(1):17.0;              
y = 2.0:h(2):17.0;
z = 0.0:h(3):12.0;
sqSz = 25;                      % for plotting

% grid struct
grid.h = h;
grid.x = x;
grid.y = y;
grid.z = z;
grid.O = [min(x), min(y), min(z)];
[X,Y,Z] = meshgrid(grid.x,grid.y,grid.z);

% DEFAULT Mallet (centered at origin)
[Phi_mallet, n_mallet] = MalletLevelSet(r,bhalf,grid,3);

% Mallet Initial Values (t=0) 
[mallet1] = MalletUp(mallet1, parabola1, t, dt, Phi_mallet, n_mallet);
[mallet2] = MalletUp(mallet2, parabola2, t, dt, Phi_mallet, n_mallet);

% Usu (static)
usu.ctr = [9.5 9.5 7]; usu.r = 6.0; usu.ht = 6.0;
[usu,cavIdx,topIdx] = UsuLevelSet(usu,grid,3);

% Colliding Grid Node indices
[Idx1] = CollisionDetect_Mallet([X(:),Y(:),Z(:)],mallet1,grid);
[Idx2] = CollisionDetect_Mallet([X(:),Y(:),Z(:)],mallet2,grid);
[Idx3] = CollisionDetect_Usu([X(:),Y(:),Z(:)],usu,grid);

% Testing: MalletUp & CollisionDetect_Mallet CollisionDetect_Usu functions
figure(1);
scatter3(X(:),Y(:),Z(:),1,'+','b'); hold on;
    % mallet
scatter3(mallet1.x_co(:,1),mallet1.x_co(:,2),mallet1.x_co(:,3),1,'ro'); hold on;
scatter3(mallet2.x_co(:,1),mallet2.x_co(:,2),mallet2.x_co(:,3),1,'co'); hold on;
    % usu
scatter3(usu.x_co(cavIdx,1),usu.x_co(cavIdx,2),usu.x_co(cavIdx,3),1,'mo'); hold on;
scatter3(usu.x_co(topIdx,1),usu.x_co(topIdx,2),usu.x_co(topIdx,3),1,'go'); hold on;
    % collided grid nodes
scatter3(X(Idx1),Y(Idx1),Z(Idx1),sqSz,'gs','filled'); hold on;
scatter3(X(Idx2),Y(Idx2),Z(Idx2),sqSz,'ms','filled'); hold on;
scatter3(X(Idx3),Y(Idx3),Z(Idx3),sqSz,'cs','filled'); hold on;
hold off;
xlabel('x'); ylabel('y'); zlabel('z');

% Track mallet centers (to see path over time)
xCent1 = mallet1.center(1); xCent2 = mallet2.center(1);
yCent1 = mallet1.center(2); yCent2 = mallet2.center(2); 
zCent1 = mallet1.center(3); zCent2 = mallet2.center(3);

while t < time(end)

    % Mallet 1
    if t >= mochiWait      
        [mallet1] = MalletUp(mallet1, parabola1, t-mochiWait, dt, Phi_mallet, n_mallet);
    end
    % Mallet 2 (wait a bit after Mallet 1)
    if t >= wait + mochiWait               
        [mallet2] = MalletUp(mallet2, parabola2, t-(mochiWait+wait), dt, Phi_mallet, n_mallet); 
    end
    
    % add new values to list of mallet centers 
    xCent1 = [xCent1; mallet1.center(1)]; xCent2 = [xCent2; mallet2.center(1)];
    yCent1 = [yCent1; mallet1.center(2)]; yCent2 = [yCent2; mallet2.center(2)]; 
    zCent1 = [zCent1; mallet1.center(3)]; zCent2 = [zCent2; mallet2.center(3)];
    
    % Mallet Speeds
    velNorm1 = norm(mallet1.v_co,2); % mallet1 
    velNorm2 = norm(mallet2.v_co,2); % mallet2 
        
    % Plot
    figure(2);
    scatter3(X(:),Y(:),Z(:),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on; % Eulerian grid
    scatter3(xCent1,yCent1,zCent1,'m'); hold on;     % mallet1 path
    scatter3(xCent2,yCent2,zCent2,'m'); hold on;     % mallet2 path
    scatter3([xCent1(end); xCent2(end)], [yCent1(end); yCent2(end)], [zCent1(end); zCent2(end)],'c','filled'); hold on; % current mallet centers
    scatter3(mallet1.x_co(:,1),mallet1.x_co(:,2),mallet1.x_co(:,3),1,'ro'); hold on; % mallet1
    scatter3(mallet2.x_co(:,1),mallet2.x_co(:,2),mallet2.x_co(:,3),1,'bo'); hold on; % mallet2
    scatter3(usu.x_co(:,1),usu.x_co(:,2),usu.x_co(:,3),1,'ko'); hold off; % usu
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([-4 23]); zlim([0 17]);
    title(['v1=',num2str(velNorm1),'   ','v2=',num2str(velNorm2),'   ','t=',num2str(t)]); 
    az = 32.7129; el = 21.1250;
    view([az,el]);
    pause(1e-10); 
    
    % update values
    k = k+1;
    t = time(k);
    
end


% =============== MALLET FUNCTIONS (KINEMATIC MALLET) =================
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
