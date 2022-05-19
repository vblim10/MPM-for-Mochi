% mochi_pound_3d_ani.m

% *** IMPORT WORKSPACE (from mochi_pound_3d.m run) ***

% VIEWS ("calibrated")
% az = 22.3778; el = 25.1439; % good angle (the usual)
% az = 0.7234; el = 90;       % top down (xy-plane)
% az = 90.0205; el = 0.0776;  % side view (yz-plane)
% az = 0.0447; el = 0.0677;   % side view (xz-plane)
az = 2.5826; el = 18.6132;  % another angle
% az = 72.5055; el = 12.0139; % angle 2
% az = 105.5813; el = 68.8487; % angle 3

% SET View (manually)
psize = 25;
sqSz = 100;
NP = size(PHX,1);
figure(1)
scatter3(PHX(:,1),PHY(:,1),PHZ(:,1),psize,'filled','c'); hold on;
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','b'); hold on;
scatter3([M1C(1,1); M2C(1,1)],[M1C(2,1); M2C(2,1)],[M1C(3,1); M2C(3,1)],'k','filled'); hold on;
scatter3(M1X(:,1),M1Y(:,1),M1Z(:,1),'ro'); hold on;
scatter3(M2X(:,1),M2Y(:,1),M2Z(:,1),'bo'); hold on;
scatter3(usu.x_co(:,1),usu.x_co(:,2),usu.x_co(:,3),1,'ko'); hold off;
xlim([-4 23]); ylim([0 19]); zlim([0 17]);
hold off;
view([az,el]);

% ANIMATION
% *** BREAK HERE *** To set VIEWING ANGLE manually ***
Jcol = JPK;
[az,el] = view;
figure(1);
set(gcf, 'color', 'w');

% gif:
% filename = 'mochi_pound_3d_ani.gif';

% mp4 video: 
filename = 'mochi_pound_3d_ani.mp4'; v = VideoWriter(filename,'MPEG-4'); open(v);

k = 1;
for n = k:length(T)
% SOLO PLOT (1 view): 
%
    scatter3(PHX(:,n),PHY(:,n),PHZ(:,n),psize*ones(NP,1),[1-Jcol(:,n),Jcol(:,n),Jcol(:,n)],'filled'); hold on; % mochi
    scatter3([M1C(1,n); M2C(1,n)],[M1C(2,n); M2C(2,n)],[M1C(3,n); M2C(3,n)],'k','filled'); hold on; % mallet centers
    scatter3(M1X(:,n),M1Y(:,n),M1Z(:,n),1,'ro'); hold on; % mallet1
    scatter3(M2X(:,n),M2Y(:,n),M2Z(:,n),1,'bo'); hold on; % mallet2
    scatter3(usu.x_co(cavIdx,1),usu.x_co(cavIdx,2),usu.x_co(cavIdx,3),1,'ko'); hold on; % usu
    xlim([-4 23]); ylim([0 19]); zlim([0 17]);
    title('\textbf{Viscoelastic: Pounded Mochi}',['$t_{k}=$',' ',num2str(T(n))],'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    hold off;
    view([az,el]);
%}
    
% DOUBLE PLOT (2 views): 
%{
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

    % Another Angle
    nexttile
    scatter3(PHX(:,n),PHY(:,n),PHZ(:,n),psize*ones(NP,1),[1-Jcol(:,n),Jcol(:,n),Jcol(:,n)],'filled'); hold on; % mochi
    scatter3([M1C(1,n); M2C(1,n)],[M1C(2,n); M2C(2,n)],[M1C(3,n); M2C(3,n)],'k','filled'); hold on; % mallet centers
    scatter3(M1X(:,n),M1Y(:,n),M1Z(:,n),1,'ro'); hold on; % mallet1
    scatter3(M2X(:,n),M2Y(:,n),M2Z(:,n),1,'bo'); hold on; % mallet2
    scatter3(usu.x_co(cavIdx,1),usu.x_co(cavIdx,2),usu.x_co(cavIdx,3),1,'ko'); hold on; % usu
    xlim([-4 23]); ylim([0 19]); zlim([0 17]);
    xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$z$','Interpreter','latex');
    az = 2.5826; el = 18.6132;  % another angle
    view([az,el]);
    
    % Angle 2
    nexttile
    scatter3(PHX(:,n),PHY(:,n),PHZ(:,n),psize*ones(NP,1),[1-Jcol(:,n),Jcol(:,n),Jcol(:,n)],'filled'); hold on; % mochi
    scatter3([M1C(1,n); M2C(1,n)],[M1C(2,n); M2C(2,n)],[M1C(3,n); M2C(3,n)],'k','filled'); hold on; % mallet centers
    scatter3(M1X(:,n),M1Y(:,n),M1Z(:,n),1,'ro'); hold on; % mallet1
    scatter3(M2X(:,n),M2Y(:,n),M2Z(:,n),1,'bo'); hold on; % mallet2
    scatter3(usu.x_co(cavIdx,1),usu.x_co(cavIdx,2),usu.x_co(cavIdx,3),1,'ko'); hold on; % usu
    xlim([-4 23]); ylim([0 19]); zlim([0 17]);
    xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex'); zlabel('$z$','Interpreter','latex');
    az = 72.5055; el = 12.0139; % angle 2
    view([az,el]);
    
    % Shared Title
    title(t,'\textbf{Viscoelastic: Pounded Mochi}',['$t_{k}=$',' ',num2str(T(n))],'Interpreter','latex');
    hold off;
%}
    
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % gif:
    %{
    if n == k
        imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append');
    end
    %}
    
    % mp4 video:
    writeVideo(v,frame);
end

% mp4 video:
close(v);

% TIMESTAMPS

% k=1;    % t=0
% k=52;   % t=0.5
% k=102;  % t=1.0
k=451;  % t=4.49

figure(2); 
Jk = JPK(:,k);
% scatter3(PHX(:,k),PHY(:,k),PHZ(:,k),psize,'filled','co'); hold on;  % original
scatter3(PHX(:,k),PHY(:,k),PHZ(:,k),psize*ones(NP,1),[1-Jk,Jk,Jk],'filled'); hold on; % mochi
scatter3([M1C(1,k); M2C(1,k)],[M1C(2,k); M2C(2,k)],[M1C(3,k); M2C(3,k)],'k','filled'); hold on; % mallet centers
scatter3(M1X(:,k),M1Y(:,k),M1Z(:,k),1,'ro'); hold on; % mallet1
scatter3(M2X(:,k),M2Y(:,k),M2Z(:,k),1,'bo'); hold on; % mallet2
scatter3(usu.x_co(cavIdx,1),usu.x_co(cavIdx,2),usu.x_co(cavIdx,3),1,'ko'); hold on; % usu
xlim([-4 23]); ylim([0 19]); zlim([0 17]);
az = 10.0953; el = 10.7466;
title('\textbf{Viscoelastic: Pounded Mochi}',['$t_{k}=$',' ',num2str(T(k))],'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
hold off;

% xlim([2 17]); ylim([3 16]); zlim([0 17]);   % v2: skinny
% az = 16.3168; el = 5.0065;  % v2: skinny

view([az,el]);
