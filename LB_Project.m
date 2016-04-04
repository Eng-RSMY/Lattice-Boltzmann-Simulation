% Lattice Bolzmann Simulation of Flow Around a Cylinder
% 
% Milestone VIII
% Bounce back around cylinder, no slip at the walls
% Grad at outlet
% Inlet is a steady flow of constant profile
% Moving Cylinder
% Grad at BB and Filling new fluid cells

clear

%% Lattice Parameters {{{1
iterations  = 10000;
pre_initialisation = 0; % # of iterations with f_eq constant in order to stabilise initial populations
Nx = 400;
Ny = 200;
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
R  = Ny/40;
Ox = Nx/4+2; Ax = 0;
Oy = Ny/2+2; Ay = 0.25;
F = 1.2; T = 570/F; % period of 1 movement cycle
o  = ( (xx-Ox).^2 + (yy-Oy).^2 ) <= R.^2; % 2D obstacle as 1s and 0s
bb = find(o); % Linear indexes of obstacle nodes
walls = (yy == 1) | (yy == Ny);
noslp = find(walls); % linear indexes of wall nodes
oi = o;
for i=1:9
    oi  = oi|circshift(o,[cx(i) cy(i)]);
end
border = find(oi-o); % linear indexes of border nodes

%% FLow Parameters {{{1
cs   = 1/sqrt(3);
U    = 0.1;
Re   = 200;
visc = 2.*U*R/Re;
tau  = visc/cs^2;
beta = 1/(2*tau+1);
alpha = 2; % initial value for the entropic over-relaxation parameter

%% Initial u, v, rho {{{1
rho = ones(Nx,Ny);
uu  = U*ones(Nx,Ny);
vv  = 0*ones(Nx,Ny);

%% Initial feq, f ,fgr {{{1
for i=1:9
  feq(i,:,:) = rho * W(i) ...
      .* (2 - sqrt(1 + 3*uu.^2)) ...
      .* (2 - sqrt(1 + 3*vv.^2)) ...
      .* ((2*uu+sqrt(1+3*uu.^2))./(1-uu)).^cx(i) ...
      .* ((2*vv+sqrt(1+3*vv.^2))./(1-vv)).^cy(i);
end
f = feq;
fgr_next = feq;

%% Initial Force Probe {{{1
Fx = 0; Fy = 0;
Ou = 0; Ov = 0; % TODO: organise
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
       tmp = rho; tmp([bb; noslp]) = nan;
       storage_rho(:,:,i) = tmp;
       tmp = uu; tmp([bb; noslp]) = nan;
       storage_uu(:,:,i) = tmp;
       tmp = vv; tmp([bb; noslp]) = nan;
       storage_vv(:,:,i) = tmp;
       % TODO: rho(bb) = nan
       storage_vprobe(i) = vv(100,50);
       storage_Fx(i) = Fx;
       storage_Fy(i) = Fy;
       storage_t(i) = t;
   end

   % Advection/Free-flight {{{2
    for i=1:9
       f(i,:,:) = circshift(f(i,:,:), [0,cx(i),cy(i)]);
    end
     
   % Boundary Conditions {{{2
    for i=1:9
        %f(i,bb) = f(opp(i),bb);       % TODO: swap would be better. but
        f(i,noslp) = f(ymr(i),noslp); % somehow MATLAB swaps automatically
    end
    for n=border'
        for i=1:9
            if o(n+cx(i)+cy(i)*Nx) % if the ith velocity leads to a solid node
                j = opp(i);
                ns = n+cx(i)+cy(i)*Nx; % that solid node's index
                f(j,n) = f(i,ns);
            end
        end
     end
    % TODO: Decide : Recalculate uu, vv, rho here?
        rhon = reshape(sum(f),Nx,Ny);
        uun = reshape((cx * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rhon;
        vvn = reshape((cy * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rhon;
    for n=border'
        % interpolate velocity derivatives
        % TODO: This doesn't work if a boundary node is trapped
        %       between two solid nodes!
        for i=2:5
            j = opp(i);
            nf = n+cx(i)+cy(i)*Nx; %guess
            if o(nf)
                nf = n+cx(j)+cy(j)*Nx; %correction
                direction = -1;
            else
                direction = 1;
            end
            if cy(i) == 0
                dudx = cx(i) * (uu(n) - uu(nf)) * direction;
                dvdx = cx(i) * (vv(n) - vv(nf)) * direction;
            end
            if cx(i) == 0
                dudy = cy(i) * (uu(n) - uu(nf)) * direction;
                dvdy = cy(i) * (vv(n) - vv(nf)) * direction;
            end
        end
        
        nDbar = 0; uutgt = 0; vvtgt = 0; rhobb = 0; rhos = 0;
        for i=1:9
            j = opp(i);
            ns = n+cx(i)+cy(i)*Nx;
            nf = n+cx(j)+cy(j)*Nx;
            if o(ns)
                nDbar = nDbar+1;
                q(j) = 1; %TODO: interpolate?
                uutgt = uutgt + (q(j)*uun(nf)+Ou)/(1+q(j));
                vvtgt = vvtgt + (q(j)*vvn(nf)+Ov)/(1+q(j));

                rhobb = rhobb + f(i,ns);
                rhos  = rhos  + 6*W(j)*rhon(n)*(cx(j)*Ou+cy(j)*Ov); % TODO: rho0 ?
            else
                rhobb = rhobb + f(j,n);
            end
        end
        uutgt = uutgt/nDbar;
        vvtgt = vvtgt/nDbar;
        rhotgt = rhobb + rhos;
        
        Pxxeq = rhotgt*cs^2 + rhotgt.*uutgt.^2;
        Pyyeq = rhotgt*cs^2 + rhotgt.*vvtgt.^2;
        Pxyeq =               rhotgt.*uutgt.*vvtgt; % = Pyxeq
        Pxx1  = -rhotgt*cs^2/2/beta * 2*dudx;
        Pyy1  = -rhotgt*cs^2/2/beta * 2*dvdy;
        Pxy1  = -rhotgt*cs^2/2/beta * (dudy + dvdx); % = Pyx1
        Pxx(n)   = Pxxeq + Pxx1;
        Pyy(n)   = Pyyeq + Pyy1;
        Pxy(n)   = Pxyeq + Pxy1; % = Pyx

        for i=1:9
            if o(n+cx(i)+cy(i)*Nx)
                j = opp(i);
                ns = n+cx(i)+cy(i)*Nx;
                nf = n+cx(j)+cy(j)*Nx;

                f(j,n) = W(j) .* ...
                        ( ...
                        + rhotgt ...
                        + rhotgt.*uutgt * cx(j) / cs^2 ...
                        + rhotgt.*vvtgt * cy(j) / cs^2 ...
                        + 1/2/cs^4 * ( (Pxx(n) - rhotgt*cs^2)*(cx(j)*cx(j) - cs^2) ...
                                      +(Pyy(n) - rhotgt*cs^2)*(cy(j)*cy(j) - cs^2) ...
                                      +2*(Pxy(n))*(cx(j)*cy(j)) ...
                                     ) ... 
                        );
 
            end
        end
     end

     % Sum up the force {{{2
     Fx = 0; Fy = 0;
     for n=border'
        for i=1:9
            if o(n+cx(i)+cy(i)*Nx) % if the ith velocity leads to a solid node
                ns = n+cx(i)+cy(i)*Nx; % that solid node's index
                Fx = Fx + cx(i)*(f(opp(i),n)+f(i,ns));
                Fy = Fy + cy(i)*(f(opp(i),n)+f(i,ns));
            end
        end
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
     Pxx  = reshape(((cx.*cx) * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
     Pyy  = reshape(((cy.*cy) * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho;
     Pxy  = reshape(((cx.*cy) * reshape(f,9,Nx*Ny)),Nx,Ny) ./ rho; % = Pyx
     fgr = fgr_next; % grad's population carried from the prev. timestep
     for i = 1:9
        fgr_next(i,:,:) = W(i) .* ...
            ( ...
            + rho ...
            + rho.*uu * cx(i) / cs^2 ...
            + rho.*vv * cy(i) / cs^2 ...
            + 1/2/cs^4 * ( (Pxx - rho*cs^2)*(cx(i)*cx(i) - cs^2) ...
                          +(Pyy - rho*cs^2)*(cy(i)*cy(i) - cs^2) ...
                          +2*(Pxy)*(cx(i)*cy(i)) ...
                         ) ... 
            );
     end
     %Cylinder
     uu(bb) = Ou;
     vv(bb) = Ov;
     % Domain
     for i = 1:9
       tmp = f(i,:,:); % save a non-collided copy of f
       
       if t >= pre_initialisation
         feq(i,:,:) = rho * W(i) ...
            .* (2 - sqrt(1 + 3*uu.^2)) ...
            .* (2 - sqrt(1 + 3*vv.^2)) ...
            .* ((2*uu+sqrt(1+3*uu.^2))./(1-uu)).^cx(i) ...
            .* ((2*vv+sqrt(1+3*vv.^2))./(1-vv)).^cy(i);
       end
       f(i,:,:) = f(i,:,:) ...
                - alpha*beta*( f(i,:,:) - feq(i,:,:) ) ;
            
       f(i,bb)     = tmp(bb);    % no collision in solid
       f(i,noslp)  = tmp(noslp); % no collision inside wall boundary
       f(i,1,2:end-1)   = feq(i,1,2:end-1);   % inlet
       if cx(i) < 0
           f(i,end,2:end-1) = fgr(i,end,2:end-1); % outlet
       end
       %f(i,bb)     = feq(i,bb); %cylinder
     end
     
   % Move Cylinder {{{2
     Ox_prev = Ox; Oy_prev = Oy;
     Ox = Nx/4+2+round(Ax*2*R*(1-cos(2*pi*t/T)));
     Oy = Ny/2+2+round(Ay*2*R*sin(2*pi*t/T));
     Ou = 2*pi/T*Ax*2*R*sin(2*pi*t/T);
     Ov = 2*pi/T*Ay*2*R*cos(2*pi*t/T); 
     o  = ( (xx-Ox).^2 + (yy-Oy).^2 ) <= R.^2; % 2D obstacle as 1s and 0s
     bb = find(o);
     oi = o;
     for i=1:9
         oi  = oi|circshift(o,[cx(i) cy(i)]);
     end
     border = find(oi-o);% Linear indexes of obstacle nodes
     % Missing populations?

   % Update Entropic Relaxation Parameter {{{2
     alpha = alpha;
     
   % Update Waitbar {{{2
   waitbar(t/iterations,w,...
       ['Simulating. - ', num2str(100*t/iterations), '% done']);
   % }}}2
end
close(w) %}}}1

%% Display {{{1

i=1;
LB_Display

% }}}1

%close all;
