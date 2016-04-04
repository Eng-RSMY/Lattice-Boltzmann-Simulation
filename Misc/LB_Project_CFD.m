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


%% Simulation {{{1
w = waitbar(0, 'Simulating.');
for t = 1:iterations
   % Storage {{{2
   if mod(t-1,storage_interval) <= 1
       if t == 1
           i = 1;
       else
           i = length(storage_t)+1 ;
       end
       clear tmp;
       tmp = rho; tmp([bb; ns]) = nan;
       storage_rho(:,:,i) = tmp;
       tmp = uu; tmp([bb; ns]) = nan;
       storage_uu(:,:,i) = tmp;
       tmp = vv; tmp([bb; ns]) = nan;
       storage_vv(:,:,i) = tmp;
       storage_t(i) = t;
   end

   % Calculations {{{2
   divrhouv = diff(rho.*uu, 1, 1)/dx + diff(rho.*vv, 1, 2)/dy;
   P = rho*cs^2  + rho.*u2 -2*nu*rho*(diff(uu,1,1)/dx+diff(vv,1,2)/dy);
   divrho
   
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
