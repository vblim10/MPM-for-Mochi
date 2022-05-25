% nH_skewColl_2d_ani.m
% Script to render simulation from nH_skewColl_2d.m as gif/mp4

% *** IMPORT WORKSPACE (from nH_skewColl_2d.m run) ***

% *** OR IMPORT DATA from .dat files (uncomment section below)***
%{
% PHX (x-coords) (NP x k+1)
    filename = 'PHX.dat'; delimiterIn = ';'; headerlinesIn = 0;
    PHX = importdata(filename,delimiterIn,headerlinesIn); 
% PHY (y-coords) (NP x k+1)
    filename = 'PHY.dat'; delimiterIn = ';'; headerlinesIn = 0;
    PHY = importdata(filename,delimiterIn,headerlinesIn); 
% JPK (NP x k+1)
    filename = 'Jpk.dat'; delimiterIn = ';'; headerlinesIn = 0;
    JPK = importdata(filename,delimiterIn,headerlinesIn);
% GNxy (grid node positions) (NG x 2)
    filename = 'GNxy.dat'; delimiterIn = ';'; headerlinesIn = 0;
    GNxy = importdata(filename,delimiterIn,headerlinesIn);
    Ei = GNxy(:,1); Ej = GNxy(:,2); 
% Energy EvT = [t,E,Q,K,El,G] (k+1 x 6)
    filename = 'Engy.dat'; delimiterIn = ';'; headerlinesIn = 0;
    EvT = importdata(filename,delimiterIn,headerlinesIn);
    T = EvT(:,1); E = EvT(:,2); El = EvT(:,3); K = EvT(:,4);
%}

% CHECK GRAPH
psize = 25; 
sqSz = 25;
NP = size(PHX,1);
figure(1)
scatter(PHX(:,1),PHY(:,1),psize,'filled','c'); hold on;
scatter(Ei(:),Ej(:),sqSz,'+','b');
hold off;

% ANIMATION
Jcol = JPK;
figure(1);
set(gcf, 'color', 'w');

% gif:
% filename = 'nH_skewColl_2d_ani.gif';

% mp4 video: 
filename = 'nH_skewColl_2d_ani.mp4'; v = VideoWriter(filename,'MPEG-4'); open(v);

for n = 1:length(T)
    % Particle color dependent on stress
    scatter(PHX(:,n),PHY(:,n),psize*ones(NP,1),[1-Jcol(:,n),Jcol(:,n),Jcol(:,n)],'filled'); hold on;
    
    % 2 color cyllinders 
%     scatter(PHX(1:NP/2,n),PHY(1:NP/2,n),psize,'filled','co'); hold on; 
%     scatter(PHX((NP/2)+1:end,n),PHY((NP/2)+1:end,n),psize,'filled','mo'); hold on; 

    % Eul. grid nodes (no color)
%     scatter(Ei(:),Ej(:),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on;

    % Eul. grid nodes (blue)
    scatter(Ei(:),Ej(:),sqSz,'+','b'); hold on; 
    
    % Title & Labels
    title('\textbf{neo-Hookean: Skew Impact of Cyllinders}',['$t_{k}=$',' ',num2str(T(n),'%.2f')],'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    hold off;
    
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
