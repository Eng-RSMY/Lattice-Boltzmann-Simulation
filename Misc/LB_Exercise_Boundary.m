clear all;
close all;

%% Iteration & Storage Parameters
iterations = 100;
N = 100;
storage_total = 10;
storage_interval = (iterations-1)/(storage_total-1);

%% LB Weights and Velocities 
W  = [ 16.0/36,  4.0/36,  4.0/36,  4.0/36,  4.0/36,  1.0/36,  1.0/36,  1.0/36,  1.0/36];
cx = [       0,       1,       0,      -1,       0,       1,      -1,      -1,       1];
cy = [       0,       0,       1,       0,      -1,       1,       1,      -1,      -1];
 
%% Physical Parameters
Re = 50;
Vmax = 0.01; %??
visc = Vmax*N/Re;
cs = 1/sqrt(3);
tau = visc/cs^2;
beta = 1/(2*tau+1) ; 

%% Boundary Nodes
[xx, yy] = meshgrid(0:N, 0:N) ;
leftboundary = (xx == 0);
rightboundary = (xx == N);
leftu = Vmax;
rightu = Vmax;

%% Initial u, v, rho
lambdax = 1; lambday = 1;
Kx = 2*pi/lambdax/N;
Ky = 2*pi/lambday/N;
K = Kx^2 + Ky^2;
Ma = Vmax/cs;
uu = - Vmax*Ky/sqrt(K) ...
    * sin(Ky*yy) .* cos(Kx*xx) ;
vv = - Vmax*Kx/sqrt(K) ...
    * sin(Kx*xx) .* cos(Ky*yy) ;
rho = 1 - Ma^2/2/K^2 * ...
    ( Ky^2 * cos(2*Kx*xx) + ...
      Kx^2 * cos(2*Ky*yy)   ); 

  uu(leftboundary) = leftu;
     uu(rightboundary) = rightu;
     vv(leftboundary) = 0;
     vv(rightboundary) = 0;
     
%% Initial f_eq
for i = 1:9
  f(:,:,i) = rho * W(i) ...
      .* (2 - sqrt(1 + 3*uu.^2)) ...
      .* (2 - sqrt(1 + 3*vv.^2)) ...
      .* ((2*uu+sqrt(1+3*uu.^2))./(1-uu)).^cx(i) ...
      .* ((2*vv+sqrt(1+3*vv.^2))./(1-vv)).^cy(i);
end


%% Simulation
w = waitbar(0, 'Simulating.');

for t = 1:iterations
    
   %storage
   if mod(t-1,storage_interval) <= 0.5
       if t == 1
           i = 1;
       else
          i = size(storage_rho,3)+1 ;
       end
       storage_rho(:,:, i) = rho;
       storage_uu(:,:, i) = uu;
       storage_vv(:,:, i) = vv;
       storage_t(i) = t;
   end
   
   %advection/free-flight
     %periodic boundary conditions
     f(:,:,2) = [f(end,:,2); f(1:end-1,:,2)];
     f(:,:,3) = [f(:,end,3), f(:,1:end-1,3)];
     f(:,:,4) = [f(2:end,:,4); f(1,:,4)    ];
     f(:,:,5) = [f(:,2:end,5), f(:,1,5)    ];
     f(:,:,6) = [f(end,end,6)  , f(end,1:end-1,6);
                 f(1:end-1,end,6), f(1:end-1,1:end-1,6)];
     f(:,:,7) = [f(end,2:end,7)  , f(end,1,7);
                 f(1:end-1,2:end,7)  , f(1:end-1,1,7)      ];
     f(:,:,8) = [f(2:end,2:end,8)  , f(2:end,1,8);
                 f(1,2:end,8)  , f(1,1,8)      ];
     f(:,:,9) = [f(2:end,end,9)  , f(2:end,1:end-1,9);
                 f(1,end,9)  , f(1,1:end-1,9)      ];
     
     %bounce-back boundary condition
%    temp = fn(1);
%    fn = [fn(2:end), fp(end)] ;
%    fp = [temp, fp(1:end-1)] ;
%    store exiting velocities
     
%     l = xx == 0;
%     r = xx == N;
%     t = yy == N;
%     b = yy == 0;
%     tl = (t+l) > 0;
%     tr = (t+r) > 0;
%     bl = (b+l) > 0;
%     br = (b+r) > 0;
     
%     a = f(:,:,2).*r;
%     b = f(:,:,3).*t;
%     c = f(:,:,4).*l;
%     d = f(:,:,5).*b;
%     e = f(:,:,6).*tr;
%     f = f(:,:,7).*tl;
%     g = f(:,:,8).*bl;
%     h = f(:,:,9).*br;
   
   
   %relaxation/collision
     rho = sum(f,3) ;
     uu = uu*0; vv = vv*0;
     for i = 1:9
         uu = uu + f(:,:,i)*cx(i);
         vv = vv + f(:,:,i)*cy(i);
     end
     uu = uu./rho; vv = vv./rho;
     
     uu(leftboundary) = leftu;
     uu(rightboundary) = rightu;
     vv(leftboundary) = 0;
     vv(rightboundary) = 0;
     
     
    
%    u = (fp - fn)./rho ;
     for i = 1:9
       feq(:,:,i) = rho * W(i) ...
            .* (2 - sqrt(1 + 3*uu.^2)) ...
            .* (2 - sqrt(1 + 3*vv.^2)) ...
            .* ((2*uu+sqrt(1+3*uu.^2))./(1-uu)).^cx(i) ...
            .* ((2*vv+sqrt(1+3*vv.^2))./(1-vv)).^cy(i);
     end

     for i = 1:9
       f(:,:,i) = f(:,:,i) ...
                - 2*beta*( f(:,:,i) - feq(:,:,i) ) ;
     end
   
    
   waitbar(t/iterations,w,...
       ['Simulating. - ', num2str(100*t/iterations), '% done']);
   
   
end
close(w)

%% Display
for i = 1:size(storage_rho,3)
    hold off
    
    scrsz = get(0,'ScreenSize');
    figure(1)
    set(1,'Name', ['t = ', num2str(storage_t(i)), ...
                   ' - Re = ', num2str(Re),...
                   ' - Press SPACE to advance'],...
          'NumberTitle', 'off')
    set(1, 'Position',[1 1 scrsz(3) scrsz(4)])
    
    subplot(1,2,1)
    surf(xx,yy,storage_rho(:,:,i))
    title('Density')
    
    subplot(1,2,2)
    quiver(xx,yy,storage_uu(:,:,i),storage_vv(:,:,i))
    title('Velocity')
    
    
    pause;
end

close all;