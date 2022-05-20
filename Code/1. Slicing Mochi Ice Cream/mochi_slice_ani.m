% mochi_slice_ani.m
% Script to render simulation from mochi_slice.m as gif/mp4

% *** IMPORT WORKSPACE (from mochi_slice.m run) ***

% VIEWS ("calibrated")
% az = 22.3778; el = 25.1439; % good angle
% az = 14.5252; el = 9.5394;  % another angle
% az = 0.7234; el = 90;       % top down (xy-plane)
% az = 90.0205; el = 0.0776;  % side view (yz-plane)
% az = 0.0447; el = 0.0677;   % side view (xz-plane)
% az = 4.8147; el = 7.2729;   % RL1 (real life view)
% az = -4.2270; el = 5.1174;  % RL2
az = -7.1795; el = 5.4121;  % RL3

% SET View (manually)
psize = 25;
NP = size(PHX,1);
figure(1)
scatter3(PHX(:,1),PHY(:,1),PHZ(:,1),psize,'filled','c'); hold on; 
scatter3(BX(:,1),BY(:,1),BZ(:,1),psize,'filled','ko'); hold on; 
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','b'); hold on;
scatter3([0],[0],[12],1,'+','r'); hold off;
xlim([0 16]); ylim([0 16]); zlim([-3 15]);
view([az,el]);

% ANIMATION
% *** BREAK HERE *** To set VIEWING ANGLE manually ***
figure(1);
[az,el] = view;
set(gcf, 'color', 'w');

% gif: 
% filename = 'mochi_slice_ani.gif';

% mp4 video:
filename = 'mochi_slice_ani.mp4'; v = VideoWriter(filename,'MPEG-4'); open(v);

for n = 1:length(T)
% SOLO PLOT (1 view):
%
    scatter3(PHX(Idx_cream,n),PHY(Idx_cream,n),PHZ(Idx_cream,n),psize,'filled','mo'); hold on;
    scatter3(PHX(Idx_mochi,n),PHY(Idx_mochi,n),PHZ(Idx_mochi,n),psize,'filled','go'); hold on;
    scatter3(BX(:,n),BY(:,n),BZ(:,n),psize,'filled','ko'); hold on;        % blade
    scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on;
    title('\textbf{Viscoelastic: Mochi Ice Cream Slice}',['$t_{k}=$',' ',num2str(T(n))],'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    xlim([0 16]); ylim([0 16]); zlim([-3 15]);
    view([az,el]);
    hold off;
%}

% DOUBLE PLOT (2 views): 
%{
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

    % RL3 view
    nexttile
    scatter3(PHX(Idx_cream,n),PHY(Idx_cream,n),PHZ(Idx_cream,n),psize,'filled','mo'); hold on; % ice cream
    scatter3(PHX(Idx_mochi,n),PHY(Idx_mochi,n),PHZ(Idx_mochi,n),psize,'filled','go'); hold on; % mochi
    scatter3(BX(:,n),BY(:,n),BZ(:,n),psize,'filled','ko'); hold on; % blade
    scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on; % grid
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    xlim([0 16]); ylim([0 16]); zlim([-3 15]);
    az = -7.1795; el = 5.4121; % RL3
    view([az,el]);    
    
    % XZ view
    nexttile
    scatter3(PHX(Idx_cream,n),PHY(Idx_cream,n),PHZ(Idx_cream,n),psize,'filled','mo'); hold on; % ice cream
    scatter3(PHX(Idx_mochi,n),PHY(Idx_mochi,n),PHZ(Idx_mochi,n),psize,'filled','go'); hold on; % mochi
    scatter3(BX(:,n),BY(:,n),BZ(:,n),psize,'filled','ko'); hold on; % blade
    scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on; % grid
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    xlim([0 16]); ylim([0 16]); zlim([-3 15]);
    az = 0.0447; el = 0.0677;   % side view (xz-plane)
    view([az,el]);  
    
    % Shared Title: 
    title(t, '\textbf{Viscoelastic: Mochi Ice Cream Slice}',['$t_{k}=$',' ',num2str(T(n))],'Interpreter','latex');
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


% TIMESTAMPS
% k=1;    % t=0
% k=30;   % t=0.4482
k=55;   % t=0.8632

figure(2)
scatter3(PHX(Idx_cream,k),PHY(Idx_cream,k),PHZ(Idx_cream,k),psize,'filled','mo'); hold on; % ice cream
scatter3(PHX(Idx_mochi,k),PHY(Idx_mochi,k),PHZ(Idx_mochi,k),psize,'filled','go'); hold on; % mochi
scatter3(BX(:,k),BY(:,k),BZ(:,k),psize,'filled','ko'); hold on; % blade
scatter3(Xg(:,1),Xg(:,2),Xg(:,3),1,'+','MarkerEdgeColor','none','MarkerFaceColor','none'); hold on; % grid
title('\textbf{Viscoelastic: Mochi Ice Cream Slice}',['$t_{k}=$',' ',num2str(T(k))],'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
xlim([0 16]); ylim([0 16]); zlim([-3 15]);

% az = 4.8147; el = 7.2729;   % RL1 (real life view)
% az = -4.2270; el = 5.1174;  % RL2
% az = -7.1795; el = 5.4121;  % RL3
% az = 22.3778; el = 25.1439; % good angle
% az = 12.7092; el = 13.5267; % another angle
% az = 90.0205; el = 0.0776;  % side view (yz-plane)
% az = 0.0447; el = 0.0677;   % side view (xz-plane)
az = 1.2360; el = 5.9727;   % RL pink

view([az,el]);



