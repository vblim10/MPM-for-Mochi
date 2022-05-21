% VE_skewColl_3d_ani.m
% Script to render simulation from VE_skewColl_3d.m as gif/mp4

% *** IMPORT WORKSPACE (from VE_skewColl_3d.m run) ***

% *** OR IMPORT DATA from .dat files (uncomment section below)***
%{
% PHX (x-coords) (NP x k+1)
    filename = 'PHX.dat'; delimiterIn = ';'; headerlinesIn = 0;
    PHX = importdata(filename,delimiterIn,headerlinesIn); 
% PHY (y-coords) (NP x k+1)
    filename = 'PHY.dat'; delimiterIn = ';'; headerlinesIn = 0;
    PHY = importdata(filename,delimiterIn,headerlinesIn); 
% PHZ (z-coords) (NP x k+1)
    filename = 'PHZ.dat'; delimiterIn = ';'; headerlinesIn = 0;
    PHZ = importdata(filename,delimiterIn,headerlinesIn); 
% JPK (NP x k+1)
    filename = 'Jpk.dat'; delimiterIn = ';'; headerlinesIn = 0;
    JPK = importdata(filename,delimiterIn,headerlinesIn);
% GNxyz (grid node positions) (NG x 3)
    filename = 'GNxyz.dat'; delimiterIn = ';'; headerlinesIn = 0;
    GNxyz = importdata(filename,delimiterIn,headerlinesIn);
    Ei = GNxyz(:,1); Ej = GNxyz(:,2); Ek = GNxyz(:,3);
% Energy EvT = [t,E,Q,K,Q_El,Q_N,G] (k+1 x 7)
    filename = 'Engy.dat'; delimiterIn = ';'; headerlinesIn = 0;
    EvT = importdata(filename,delimiterIn,headerlinesIn);
    T = EvT(:,1); E = EvT(:,2); El = EvT(:,3); K = EvT(:,4);
%}

% VIEWS ("calibrated")
% az = 22.3778; el = 25.1439; % good angle
az = 8.3243; el = 25.1439;  % good angle 2
% az = 0.7234; el = 90;       % top down (xy-plane)
% az = 90.0205; el = 0.0776;  % side view (yz-plane)
% az = 0.0447; el = 0.0677;   % side view (xz-plane)
view([az,el]);

% SET View (manually)
psize = 25; 
NP = size(PHX,1);
figure(1)
scatter3(PHX(:,1),PHY(:,1),PHZ(:,1),psize,'filled','c'); hold on;
scatter3(Ei(:),Ej(:),Ek(:),1,'+','b');
hold off;
view([az,el]);

% ANIMATION
% *** BREAK HERE *** To set VIEWING ANGLE manually ***
Jcol = JPK;
[az,el] = view;
figure(1);
set(gcf, 'color', 'w');

% gif:
% filename = 'VE_skewColl_3d_ani.gif';

% mp4 video: 
filename = 'VE_skewColl_3d_ani.mp4'; v = VideoWriter(filename,'MPEG-4'); open(v);

for n = 1:length(T)
    % Particle color dependent on stress
%     scatter3(PHX(:,n),PHY(:,n),PHZ(:,n),psize*ones(NP,1),[1-Jcol(:,n),Jcol(:,n),Jcol(:,n)],'filled'); hold on;
    
    % 2 color balls 
    scatter3(PHX(1:NP/2,n),PHY(1:NP/2,n),PHZ(1:NP/2,n),psize,'filled','co'); hold on; 
    scatter3(PHX((NP/2)+1:end,n),PHY((NP/2)+1:end,n),PHZ((NP/2)+1:end,n),psize,'filled','mo'); hold on; 

    % Eul. grid
    scatter3(Ei(:),Ej(:),Ek(:),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold off;
    
    % Title & Labels
    title('\textbf{Viscoelastic: Skew Impact of Balls}',['$t_{k}=$',' ',num2str(T(n),'%.2f')],'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    
    view([az,el]);
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % gif:
    %{
    if n == 1
        imwrite(imind,cm,filename,'gif','DelayTime',0.00001, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.00001,'WriteMode','append');
    end
    %}
    
    % mp4 video:
    writeVideo(v,frame);
end

% mp4 video:
close(v);
