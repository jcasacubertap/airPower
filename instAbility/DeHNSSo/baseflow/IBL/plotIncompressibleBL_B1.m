%% Inputs

%Define 4 main x-stations (array position) to evaluate BL profiles
plt.st1 = 1;
plt.st2 = round(length(BL.x(1,1:end))/3);
plt.st3 = round(length(BL.x(1,1:end))/3*2);
plt.st4 = length(BL.x(1,1:end));

%Define 6 alternative x-stations (array position) to evaluate BL profiles
plt.st1Alt = 1;
plt.st2Alt = round(length(BL.x(1,1:end))/5);
plt.st3Alt = round(length(BL.x(1,1:end))/5*2);
plt.st4Alt = round(length(BL.x(1,1:end))/5*3);
plt.st5Alt = round(length(BL.x(1,1:end))/5*4);
plt.st6Alt = length(BL.x(1,1:end));

%Change to smoother plotting capabilities of matlab
set(0, 'DefaultFigureRenderer', 'painters');

%% Plot main fields 
%% u-velocity

hFig = figure(1);
set(hFig,'Position',[200,600,1400,700]);

%Velocity full colormap
subplot(2,4,[1 2])
hold on
pcolor(BL.x,BL.y,BL.u); shading interp; colormap('jet'); colorbar;
xlabel('$x \; (\textrm{m})$','interpreter','latex');
ylabel('$y \; (\textrm{m})$','interpreter','latex');
title('$u \; (\textrm{m/s})$','interpreter','latex');
caxis([min(BL.u(:)) max(BL.u(:))]);
axis([BL.x(1) BL.x(end) BL.y(end) 2*BL.y_i]);

%Profiles at fixed x-position
subplot(2,4,3)
hold on; grid on;
plot(BL.u(:,plt.st1Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$u$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st1Alt))],'interpreter','latex')
axis([min(BL.u(:,plt.st1Alt)) BL.u(1,plt.st1Alt)*1.2 0 2*BL.y_i]);

subplot(2,4,4)
hold on; grid on;
plot(BL.u(:,plt.st2Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$u$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st2Alt))],'interpreter','latex')
axis([min(BL.u(:,plt.st2Alt)) BL.u(1,plt.st2Alt)*1.2 0 2*BL.y_i]);

subplot(2,4,5)
hold on; grid on;
plot(BL.u(:,plt.st3Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$u$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st3Alt))],'interpreter','latex')
axis([min(BL.u(:,plt.st3Alt)) BL.u(1,plt.st3Alt)*1.2 0 2*BL.y_i]);

subplot(2,4,6)
hold on; grid on;
plot(BL.u(:,plt.st4Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$u$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st4Alt))],'interpreter','latex')
axis([min(BL.u(:,plt.st4Alt)) BL.u(1,plt.st4Alt)*1.2 0 2*BL.y_i]);

subplot(2,4,7)
hold on; grid on;
plot(BL.u(:,plt.st5Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$u$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st5Alt))],'interpreter','latex')
axis([min(BL.u(:,plt.st5Alt)) BL.u(1,plt.st5Alt)*1.2 0 2*BL.y_i]);

subplot(2,4,8)
hold on; grid on;
plot(BL.u(:,plt.st6Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$u$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st6Alt))],'interpreter','latex')
axis([min(BL.u(:,plt.st6Alt)) BL.u(1,plt.st6Alt)*1.2 0 2*BL.y_i]);

%% v-velocity

hFig = figure(2);
set(hFig,'Position',[200,600,1400,700]);

%Velocity full colormap
subplot(2,4,[1 2])
hold on
pcolor(BL.x,BL.y,BL.v); shading interp; colormap('jet'); colorbar;
xlabel('$x \; (\textrm{m})$','interpreter','latex');
ylabel('$y \; (\textrm{m})$','interpreter','latex');
title('$v \; (\textrm{m/s})$','interpreter','latex');
caxis([min(BL.v(:)) max(BL.v(:))]);
axis([BL.x(1) BL.x(end) BL.y(end) 2*BL.y_i]);

%Profiles at fixed x-position
subplot(2,4,3)
hold on; grid on;
plot(BL.v(:,plt.st1Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$v$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st1Alt))],'interpreter','latex')
axis([min(BL.v(:,plt.st1Alt))*1.2 max(BL.v(:,plt.st1Alt))*1.2 0 2*BL.y_i]);

subplot(2,4,4)
hold on; grid on;
plot(BL.v(:,plt.st2Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$v$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st2Alt))],'interpreter','latex')
axis([min(BL.v(:,plt.st2Alt))*1.2 max(BL.v(:,plt.st2Alt))*1.2 0 2*BL.y_i]);

subplot(2,4,5)
hold on; grid on;
plot(BL.v(:,plt.st3Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$v$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st3Alt))],'interpreter','latex')
axis([min(BL.v(:,plt.st3Alt))*1.2 max(BL.v(:,plt.st3Alt))*1.2 0 2*BL.y_i]);

subplot(2,4,6)
hold on; grid on;
plot(BL.v(:,plt.st4Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$v$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st4Alt))],'interpreter','latex')
axis([min(BL.v(:,plt.st4Alt))*1.2 max(BL.v(:,plt.st4Alt))*1.2 0 2*BL.y_i]);

subplot(2,4,7)
hold on; grid on;
plot(BL.v(:,plt.st5Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$v$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st5Alt))],'interpreter','latex')
axis([min(BL.v(:,plt.st5Alt))*1.2 max(BL.v(:,plt.st5Alt))*1.2 0 2*BL.y_i]);

subplot(2,4,8)
hold on; grid on;
plot(BL.v(:,plt.st6Alt),BL.y,'o-','Color','k','LineWidth',1.2,...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
xlabel('$v$','interpreter','latex')
ylabel('$y$','interpreter','latex')
title(['$x = \;$' num2str(BL.x(1,plt.st6Alt))],'interpreter','latex')
axis([min(BL.v(:,plt.st6Alt))*1.2 max(BL.v(:,plt.st6Alt))*1.2 0 2*BL.y_i]);

%% w-velocity

if BL.We ~= 0

    hFig = figure(3);
    set(hFig,'Position',[200,600,1400,700]);
    
    %Velocity full colormap
    subplot(2,4,[1 2])
    hold on
    pcolor(BL.x,BL.y,abs(BL.w)); shading interp; colormap('jet'); colorbar;
    xlabel('$x \; (\textrm{m})$','interpreter','latex');
    ylabel('$y \; (\textrm{m})$','interpreter','latex');
    title('$|w| \; (\textrm{m/s})$','interpreter','latex');
    caxis([min(BL.w(:)) max(BL.w(:))]);
    axis([BL.x(1) BL.x(end) BL.y(end) 2*BL.y_i]);
    
    %Profiles at fixed x-position
    subplot(2,4,3)
    hold on; grid on;
    plot(abs(BL.w(:,plt.st1Alt)),BL.y,'o-','Color','k','LineWidth',1.2,...
            'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    xlabel('$|w|$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$x = \;$' num2str(BL.x(1,plt.st1Alt))],'interpreter','latex')
    axis([min(abs(BL.w(:,plt.st1Alt))) abs(BL.w(1,plt.st1Alt))*1.2 0 2*BL.y_i]);
    
    subplot(2,4,4)
    hold on; grid on;
    plot(abs(BL.w(:,plt.st2Alt)),BL.y,'o-','Color','k','LineWidth',1.2,...
            'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    xlabel('$|w|$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$x = \;$' num2str(BL.x(1,plt.st2Alt))],'interpreter','latex')
    axis([min(abs(BL.w(:,plt.st2Alt))) abs(BL.w(1,plt.st2Alt))*1.2 0 2*BL.y_i]);
    
    subplot(2,4,5)
    hold on; grid on;
    plot(abs(BL.w(:,plt.st3Alt)),BL.y,'o-','Color','k','LineWidth',1.2,...
            'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    xlabel('$|w|$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$x = \;$' num2str(BL.x(1,plt.st3Alt))],'interpreter','latex')
    axis([min(abs(BL.w(:,plt.st3Alt))) abs(BL.w(1,plt.st3Alt))*1.2 0 2*BL.y_i]);
    
    subplot(2,4,6)
    hold on; grid on;
    plot(abs(BL.w(:,plt.st4Alt)),BL.y,'o-','Color','k','LineWidth',1.2,...
            'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    xlabel('$|w|$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$x = \;$' num2str(BL.x(1,plt.st4Alt))],'interpreter','latex')
    axis([min(abs(BL.w(:,plt.st4Alt))) abs(BL.w(1,plt.st4Alt))*1.2 0 2*BL.y_i]);
    
    subplot(2,4,7)
    hold on; grid on;
    plot(abs(BL.w(:,plt.st5Alt)),BL.y,'o-','Color','k','LineWidth',1.2,...
            'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    xlabel('$|w|$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$x = \;$' num2str(BL.x(1,plt.st5Alt))],'interpreter','latex')
    axis([min(abs(BL.w(:,plt.st5Alt))) abs(BL.w(1,plt.st5Alt))*1.2 0 2*BL.y_i]);
    
    subplot(2,4,8)
    hold on; grid on;
    plot(abs(BL.w(:,plt.st6Alt)),BL.y,'o-','Color','k','LineWidth',1.2,...
            'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
    xlabel('$|w|$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['$x = \;$' num2str(BL.x(1,plt.st6Alt))],'interpreter','latex')
    axis([min(abs(BL.w(:,plt.st6Alt))) abs(BL.w(1,plt.st6Alt))*1.2 0 2*BL.y_i]);

end

%% Plot inviscid streamlines

%Define conditions for streamline integration
aux.nz            = 12; %number of streamlines
aux.ylimStream    = 2;  %length in spanwise direction
[aux.xsp,aux.zsp] = meshgrid(BL.x(1,:),linspace(aux.ylimStream,-aux.ylimStream,aux.nz));
aux.usp           = repmat(BL.u(1,:),aux.nz,1);
aux.wsp           = repmat(BL.w(1,:),aux.nz,1);

%Compute and plot streamlines
hFig = figure(4);
set(hFig,'Position',[200,600,1400,600]);
sgtitle(['Inviscid streamlines'],'interpreter','latex');

hlines1 = streamline(aux.xsp,aux.zsp,aux.usp,aux.wsp,aux.xsp(:,1),aux.zsp(:,1),[1,10000]);
set(hlines1,'LineWidth',1.2,'Color','b'); grid on; set(gca,'ydir','reverse');

%Plot surface for visalization purposes
line([aux.xsp(1,1) aux.xsp(1,end)],[-aux.ylimStream -aux.ylimStream],'LineWidth',1.5,'Color','k','LineStyle','-')
hold on
line([aux.xsp(1,1) aux.xsp(1,end)],[aux.ylimStream aux.ylimStream],'LineWidth',1.5,'Color','k','LineStyle','-')
line([aux.xsp(1,1) aux.xsp(1,1)],[-aux.ylimStream -aux.ylimStream],'LineWidth',1.5,'Color','k','LineStyle','-')
p = patch([aux.xsp(1,1) aux.xsp(1,1) aux.xsp(1,end) aux.xsp(1,end)],[-aux.ylimStream aux.ylimStream aux.ylimStream -aux.ylimStream],'k');
set(p,'FaceAlpha',0.2); xlabel('$x \; (\textrm{m})$','interpreter','latex'); ylabel('$z \; (\textrm{m})$','interpreter','latex');

clear aux

%% Plot edge, displacement, and momentum thickness metrics

hFig = figure(5);
sgtitle(['Boundary Layer properties'],'interpreter','latex');
set(hFig,'Position',[200,600,800,600]);

subplot(3,1,1)
plot(BL.x,BL.edge,'k-','LineWidth',1.6); grid on;
xlabel('$x \; (\textrm{m})$','interpreter','latex');
ylabel('$\delta_{99} \; (\textrm{m})$','interpreter','latex');

subplot(3,1,2)
plot(BL.x,BL.disp,'b-','LineWidth',1.6); grid on;
xlabel('$x \; (\textrm{m})$','interpreter','latex');
ylabel('$\delta^{\star} \; (\textrm{m})$','interpreter','latex');

subplot(3,1,3)
plot(BL.x,BL.mom,'m-','LineWidth',1.6); grid on;
xlabel('$x \; (\textrm{m})$','interpreter','latex');
ylabel('$\theta \; (\textrm{m})$','interpreter','latex');

