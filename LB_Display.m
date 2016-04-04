maxrho = max(max(max(storage_rho)));
minrho = min(min(min(storage_rho)));
if minrho >= maxrho
    maxrho = 1.1*minrho;
end
storage_u2 = sqrt(storage_uu.^2+storage_vv.^2);
Cd = storage_Fx/0.5/U^2/ceil(2*R);



scrsz = get(0,'ScreenSize');
    figure(1)
    set(1, 'Position',[1 1 scrsz(3) scrsz(4)])
   
    pause_btn = uicontrol('Style', 'pushbutton', 'String', 'Play/Pause',...
        'Position', [20 20 100 40],...
        'Callback', 'if running==0, LB_Display, else, running=0;,end'); 
    
    subplot(2,1,1)
        cyl_plot = cylinder3D([Ox Oy],[minrho maxrho],R,10);
        hold on
        rho_plot = surf(xx,yy,storage_rho(:,:,i),...
            'edgecolor','none');
        hold off
        title('Density')
        zlim([minrho maxrho])
        view(0,0)
        
        subplot(2,1,2)
        u2_plot  = surf(xx,yy,...
            storage_u2(:,:,i),...
            'edgecolor', 'none');
        axis([0 Nx -0.001 Ny])    
        view(0,90)
        hold on
        n = Ny/5;
        uv_plot = quiver3(xx(1:n:end,1:n:end),...
            yy(1:n:end,1:n:end),...
            ones(size(xx(1:n:end,1:n:end))),...
            storage_uu(1:n:end,1:n:end,i),...
            storage_vv(1:n:end,1:n:end,i),...
            zeros(size(xx(1:n:end,1:n:end))));
        title('Velocity')
        hold off

        set(rho_plot, 'ZDataSource', 'storage_rho(:,:,i)')
        set(u2_plot,  'ZDataSource', 'storage_u2(:,:,i)')
        set(uv_plot,  'UDataSource', 'storage_uu(1:n:end,1:n:end,i)')
        set(uv_plot,  'VDataSource', 'storage_vv(1:n:end,1:n:end,i)')
% Animation
if exist('i')==0
  i=0;
end
running = 1;
while running
    figure(1)
    set(1,'Name', ['t = ', num2str(storage_t(i)), ...
                   ' - Re = ', num2str(Re),...
                   ' - Press SPACE to advance'],...
          'NumberTitle', 'off')
    
    refreshdata(rho_plot)
    refreshdata(u2_plot)
    refreshdata(uv_plot)

    pause(0.02);
    
    i = i+1;
    if i > size(storage_rho,3)
        i = 1;
    end
end
