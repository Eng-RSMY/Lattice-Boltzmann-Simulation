iterations = 3000;
N = 800;
storage_total = 2;
storage_interval = iterations/(storage_total-1);
clear storage;

%c = 1
%dt = 1
%dx = 1


cs = 1/sqrt(3);
u = 0.1; %is u constant?
D = 5*10^-5 ;
beta = 1/(1+6*D) ; % More on this please

x = 0:N ;
%rho = 1 + 0.5 * exp(-5000*(x/N - 1/4).^2 ) ; %steep gaussian
delta = 0.01;
rho = 1 + 0.5 * (1 - tanh((x/N - 2/10)/delta));
%rho = 1 + 0.5 * ( (N/4 < x).* (x < N/2) ) ;

fz = 2*rho/3 .* (2 - sqrt(1 + u.^2/cs^2)) ; 
fn =   rho/3 .* ((-u-cs^2)/2/cs^2 + sqrt(1 + u.^2/cs^2)) ;
fp =   rho/3 .* ((u-cs^2)/2/cs^2 + sqrt(1 + u.^2/cs^2)) ;
storage(1,:) = rho;

for t = 1:iterations
   %advection/free-flight
     %periodic boundary conditions
     fn = [fn(2:end), fn(1)] ;
     fp = [fp(end), fp(1:end-1)] ;
   
     %bounce-back boundary condition
%    temp = fn(1);
%    fn = [fn(2:end), fp(end)] ;
%    fp = [temp, fp(1:end-1)] ;
   
   
   %relaxation/collision
     rho = fz + fp + fn ;
%    u = (fp - fn)./rho ;
     fzeq = 2*rho/3 .* (2 - sqrt(1 + u.^2/cs^2)) ;  
     fneq =   rho/3 .* ((-u-cs^2)/2/cs^2 + sqrt(1 + u.^2/cs^2)) ;
     fpeq =   rho/3 .* ((u-cs^2)/2/cs^2 + sqrt(1 + u.^2/cs^2)) ;
    
     fz = fz - 2*beta*( fz - fzeq ) ;
     fn = fn - 2*beta*( fn - fneq ) ;
     fp = fp - 2*beta*( fp - fpeq ) ;
    
    
    if mod(t,storage_interval) == 0
        i = 1+(t/storage_interval) ;
        storage(i,:) = rho;
    end
    
end
    
for i = 1:1+(iterations/storage_interval)
    hold all;
    plot(storage(i,:))
end