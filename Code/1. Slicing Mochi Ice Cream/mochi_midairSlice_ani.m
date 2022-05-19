% mochi_midairSlice_ani.m

% *** IMPORT WORKSPACE (from mochi_midairSlice.m run) ***

% VIEWS ("calibrated")
% az = 22.3778; el = 25.1439; % v1
az = 5.6770; el = 9.7807;   % v2
% az = 22.3778; el = 25.1439; % good angle
% az = 14.5252; el = 9.5394;  % another angle
% az = 0.7234; el = 90;       % top down (xy-plane)
% az = 90.0205; el = 0.0776;  % side view (yz-plane)
% az = 0.0447; el = 0.0677;   % side view (xz-plane)
% az = 4.8147; el = 7.2729;   % RL1 (real life view)
% az = -4.2270; el = 5.1174;  % RL2
% az = -7.1795; el = 5.4121;  % RL3

% SET View (manually)
psize = 25;
NP = size(PHX,1);
figure(1)
scatter3(PHX(:,1),PHY(:,1),PHZ(:,1),psize,'filled','c'); hold on;
scatter3(B1X(:,1),B1Y(:,1),B1Z(:,1),psize,'filled','ko'); hold on; 
scatter3(B2X(:,1),B2Y(:,1),B2Z(:,1),psize,'filled','ko'); hold on; 
scatter3(B3X(:,1),B3Y(:,1),B3Z(:,1),psize,'filled','ko'); hold on; 
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','b'); hold on;
xlim([0 16]); ylim([0 16]); zlim([0 14]); hold off;
view([az,el]);

% ANIMATION
% *** BREAK HERE *** To set VIEWING ANGLE manually ***
figure(1);
[az,el] = view;
set(gcf, 'color', 'w');

% gif:
% filename = 'mochi_midairSlice_ani.gif';

% mp4 video: 
filename = 'mochi_midairSlice_ani.mp4'; v = VideoWriter(filename,'MPEG-4'); open(v);

for n = 1:length(T)
% SOLO PLOT (1 view):
%{
    scatter3(PHX(Idx_cream,n),PHY(Idx_cream,n),PHZ(Idx_cream,n),psize,'filled','mo'); hold on; % ice cream
    scatter3(PHX(Idx_mochi,n),PHY(Idx_mochi,n),PHZ(Idx_mochi,n),psize,'filled','go'); hold on; % mochi
    scatter3(B1X(:,n),B1Y(:,n),B1Z(:,n),psize,'filled','ko'); hold on; % blade 1
    scatter3(B2X(:,n),B2Y(:,n),B2Z(:,n),psize,'filled','ko'); hold on; % blade 2
    scatter3(B3X(:,n),B3Y(:,n),B3Z(:,n),psize,'filled','ko'); hold on; % blade 3
    xlim([0 16]); ylim([0 16]); zlim([0 14]); 
    title('\textbf{Viscoelastic: Midair Mochi Ice Cream Slice}',['$t_{k}=$',' ',num2str(T(n))],'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    hold off;
    view([az,el]);
%}

% DOUBLE PLOT (2 views): 
%
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

    % View 2
    nexttile
    scatter3(PHX(Idx_cream,n),PHY(Idx_cream,n),PHZ(Idx_cream,n),psize,'filled','mo'); hold on; % ice cream
    scatter3(PHX(Idx_mochi,n),PHY(Idx_mochi,n),PHZ(Idx_mochi,n),psize,'filled','go'); hold on; % mochi
    scatter3(B1X(:,n),B1Y(:,n),B1Z(:,n),psize,'filled','ko'); hold on; % blade 1
    scatter3(B2X(:,n),B2Y(:,n),B2Z(:,n),psize,'filled','ko'); hold on; % blade 2
    scatter3(B3X(:,n),B3Y(:,n),B3Z(:,n),psize,'filled','ko'); hold on; % blade 3
    xlim([0 16]); ylim([0 16]); zlim([0 14]); 
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    az = 5.6770; el = 9.7807; % v2
    view([az,el]);
    
    % xz-plane
    nexttile
    scatter3(PHX(Idx_cream,n),PHY(Idx_cream,n),PHZ(Idx_cream,n),psize,'filled','mo'); hold on; % ice cream
    scatter3(PHX(Idx_mochi,n),PHY(Idx_mochi,n),PHZ(Idx_mochi,n),psize,'filled','go'); hold on; % mochi
    scatter3(B1X(:,n),B1Y(:,n),B1Z(:,n),psize,'filled','ko'); hold on; % blade 1
    scatter3(B2X(:,n),B2Y(:,n),B2Z(:,n),psize,'filled','ko'); hold on; % blade 2
    scatter3(B3X(:,n),B3Y(:,n),B3Z(:,n),psize,'filled','ko'); hold on; % blade 3
    xlim([0 16]); ylim([0 16]); zlim([0 14]); 
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    az = 0.0447; el = 0.0677;   % side view (xz-plane)
    view([az,el]);
    
    % Shared Title
    title(t,'\textbf{Viscoelastic: Midair Mochi Ice Cream Slice}',['$t_{k}=$',' ',num2str(T(n))],'Interpreter','latex');
    hold off;
%}

    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % gif: 
    %{
    if n == 1
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

% Timestamps
% k=1;    % t=0
% k=30;   % t=0.4482
% k=55;   % t=0.8632
k=100;  % t=1.6102
% k=140;  % t=2.2742

figure(2)
scatter3(PHX(Idx_cream,k),PHY(Idx_cream,k),PHZ(Idx_cream,k),psize,'filled','mo'); hold on; % ice cream
scatter3(PHX(Idx_mochi,k),PHY(Idx_mochi,k),PHZ(Idx_mochi,k),psize,'filled','go'); hold on; % mochi
scatter3(B1X(:,k),B1Y(:,k),B1Z(:,k),psize,'filled','ko'); hold on; % blade 1
scatter3(B2X(:,k),B2Y(:,k),B2Z(:,k),psize,'filled','ko'); hold on; % blade 2
scatter3(B3X(:,k),B3Y(:,k),B3Z(:,k),psize,'filled','ko'); hold on; % blade 3
title('\textbf{Viscoelastic: Midair Mochi Ice Cream Slice}',['$t_{k}=$',' ',num2str(T(k))],'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
xlim([0 16]); ylim([0 16]); zlim([2 16]); 

% az = 22.3778; el = 25.1439; % good angle
% az = 12.7092; el = 13.5267; % another angle
% az = 90.0205; el = 0.0776;  % side view (yz-plane)
az = 0.0447; el = 0.0677;   % side view (xz-plane)

view([az,el]);


