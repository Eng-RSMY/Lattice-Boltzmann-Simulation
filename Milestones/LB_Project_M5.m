% Lattice Bolzmann Simulation of Flow Around a Cylinder
% 
% Bounce back around cylinder, no slip at the walls
% Grad at outlet
% Inlet is a steady flow of constant profile

clear

%% Lattice Parameters {{{1
iterations  = 400;
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
storage_total = 400; % Amount of timepoints at which values are recorded
storage_interval = (iterations-1)/(storage_total-1);

%% Obstacle & Wall Parameters {{{1
R  = Ny/30;
Ox = Nx/4+2;
Oy = Ny/2+2;
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
Re   = 1;
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

%% Initial Force Probe {{{1
Fx = 0; Fy = 0;
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
     
   % Boundary Conditions {{{2
    for i=1:9
        f(i,bb) = f(opp(i),bb); % TODO: swap is better.
        f(i,noslp) = f(ymr(i),noslp);
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
       f(i,noslp)  = tmp(noslp);
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

i=1;
LB_Display

% }}}1

close all;
