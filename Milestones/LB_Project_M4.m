% Lattice Bolzmann Simulation of Flow Around a Cylinder
% 
% Bounce back around cylinder, no slip at the walls
% Grad at outlet
% Inlet is a steady flow of constant profile

clear

%% Lattice Parameters {{{1
iterations  = 10000;
pre_initialisation = 0; % # of iterations with f_eq constant in order to stabilise initial populations
Nx = 200;
Ny = 100;
x  = [1:Nx];  % Cells in x-direction
y  = [1:Ny];  % Cells in y-direction
[xx,yy] = meshgrid(x,y);
xx = xx'; yy = yy';
W  = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36];
% i     1,   2,  3,  4,  5,    6,   7,   8,   9 
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1];
% opp contains the indices of opposite velocities
% i.e: cx(opp(i)) = -cx(i),    cy(opp(i)) = -cy(i)
opp= [  1,   4,  5,  2,  3,    8,   9,   6,   7];
% ymr contains the indices of y-mirrored velocities
ymr= [  1,   2,  5,  4,  3,    9,   8,   7,   6];

%% Storage Parameters {{{1
storage_total = 1000; % Amount of timepoints at which values are recorded
storage_interval = (iterations-1)/(storage_total-1);

%% Obstacle & Wall Parameters {{{1
R  = Ny/30;
Ox = Nx/4+2;
Oy = Ny/2+2;
o  = ( (xx-Ox).^2 + (yy-Oy).^2 ) <= R.^2; % 2D obstacle as 1s and 0s
bb = find(o); % Linear indexes of obstacle nodes
walls = (yy == 1) | (yy == Ny);
ns = find(walls); % linear indexes of wall nodes

%% FLow Parameters {{{1
cs   = 1/sqrt(3);
U    = 0.1;
Re   = 50;
visc = 2.*U*R/Re;
tau  = visc/cs^2;
beta = 1/(2*tau+1);
alpha = 2; % initial value for the entropic over-relaxation parameter

%% Initial u, v, rho {{{1
rho = ones(Nx,Ny);
uu  = U*ones(Nx,Ny);
vv  = 0*ones(Nx,Ny);

%% Initial f_eq, f {{{1
for i=1:9
  feq(i,:,:) = rho * W(i) ...
      .* (2 - sqrt(1 + 3*uu.^2)) ...
      .* (2 - sqrt(1 + 3*vv.^2)) ...
      .* ((2*uu+sqrt(1+3*uu.^2))./(1-uu)).^cx(i) ...
      .* ((2*vv+sqrt(1+3*vv.^2))./(1-vv)).^cy(i);
end
f = feq;
% }}}1

%% Simulation {{{1
w = waitbar(0, 'Simulating.');
for t = 1:iterations
   % Storage {{{2
   if mod(t-1,storage_interval) <= 1
       if t == 1
           i = 1;
       else
          i = size(storage_rho,3)+1 ;
       end
       clear tmp;
       tmp = rho; tmp([bb; ns]) = nan;
       storage_rho(:,:,i) = tmp;
       tmp = uu; tmp([bb; ns]) = nan;
       storage_uu(:,:,i) = tmp;
       tmp = vv; tmp([bb; ns]) = nan;
       storage_vv(:,:,i) = tmp;
       % TODO: rho(bb) = nan
       storage_t(i) = t;
   end

   % Advection/Free-flight {{{2
    for i=1:9
       f(i,:,:) = circshift(f(i,:,:), [0,cx(i),cy(i)]);
    end

   % Boundary Conditions {{{2
    for i=1:9
        f(i,bb) = f(opp(i),bb); % TODO: swap is better.
        f(i,ns) = f(ymr(i),ns);
    end

   % Relaxation/Collision {{{2
     rho = reshape(sum(f),Nx,Ny);
     uu = reshape((cx * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
     vv = reshape((cy * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
     % Inlet
     uu(1,2:end-1) = U;
     vv(1,2:end-1) = 0;
     rho(1,2:end-1) = 1; % rho before advec -> constant
     % Outlet
     % Domain
     for i = 1:9
       if t >= pre_initialisation
         feq(i,:,:) = rho * W(i) ...
            .* (2 - sqrt(1 + 3*uu.^2)) ...
            .* (2 - sqrt(1 + 3*vv.^2)) ...
            .* ((2*uu+sqrt(1+3*uu.^2))./(1-uu)).^cx(i) ...
            .* ((2*vv+sqrt(1+3*vv.^2))./(1-vv)).^cy(i);
       end
       tmp = f(i,:,:);
       f(i,:,:) = f(i,:,:) ...
                - alpha*beta*( f(i,:,:) - feq(i,:,:) ) ;
       f(i,bb)  = tmp(bb);
       f(i,ns)  = tmp(ns);
       f(i,1,2:end-1) = feq(i,1,2:end-1);
     end


   % Update Entropic Relaxation Parameter {{{2
     alpha = alpha;

   % Update Waitbar {{{2
   waitbar(t/iterations,w,...
       ['Simulating. - ', num2str(100*t/iterations), '% done']);
   % }}}2
end
close(w) %}}}1

%% Display {{{1
maxrho = max(max(max(storage_rho)));
minrho = min(min(min(storage_rho)));
storage_u2 = sqrt(storage_uu.^2+storage_vv.^2);
for i = 1:size(storage_rho,3)
    hold off
    
    scrsz = get(0,'ScreenSize');
    figure(1)
    set(1,'Name', ['t = ', num2str(storage_t(i)), ...
                   ' - Re = ', num2str(Re),...
                   ' - Press SPACE to advance'],...
          'NumberTitle', 'off')
    set(1, 'Position',[1 1 scrsz(3) scrsz(4)])
    
    if i==1
        subplot(2,1,1)
        rho_plot = surf(xx,yy,storage_rho(:,:,i),...
            'edgecolor','none');
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
    end
    refreshdata(rho_plot)
    refreshdata(u2_plot)
    refreshdata(uv_plot)

    pause(0.02);
end
% }}}1

close all;
