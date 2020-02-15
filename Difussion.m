%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Matlab code for Poisson equation and diffusion equation      %
%                  Zhenyang Yuan                %
%                   12/02/2020                  %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all

%% initialization
L = 1;%%length of the rod
n = 10;%number of elements
t_total = 0.5;

alpha = 1;

%%
x = linspace(-L/2,L/2,n+1);


%%
x(2:length(x)-1) = x(2:length(x)-1) + (-1)^round(rand(1))*rand(1,length(x)-2)*0.5*(x(2)-x(1));    %%rand mesh size
dx(1:length(x)-1) = x(2:length(x)) - x(1:length(x)-1) ;

cfl = 0.17;      %%fake CFL   test for optime one, this one is ideal for n = 10;
dt = (min(dx)/2+max(dx)/2)^2*cfl; 
ndt = t_total/dt;

%%construct matrix L_e and f_e
L_e = cell(1,n);
for i = 1 : n
        L_e(:,i) = {[1,-1;-1,1]/dx(i)};%
end

f_e = cell(1,n);
for i = 1 : n
        f_e(:,i) = {[ 1*100/dx(i)*(exp(x(i+1)) - (dx(i)+1)*exp(x(i)));...
                        1*100/dx(i)*((dx(i)-1)*exp(x(i+1))+exp(x(i)))]};
end


M_e = cell(1,n);
for i = 1 : n
        M_e(:,i) = {[1/3,1/6;1/6,1/3]*dx(i)};%
end

%%matrix A
A = sparse(2*n,n+1);
A_eye = [1;1];
B_eye = eye(n);
A_t = kron(B_eye,A_eye);
A(2:2*n,2:n+1) = A_t(1:2*n-1,1:n);
A(1,1) = 1;

%%extend matrix L and f and M
LL = sparse(2*n,2*n);
FF = sparse(2*n,1);
MM = sparse(2*n,2*n);
for i = 1 : n
    LL(2*i-1:2*i,2*i-1:2*i) = cell2mat(L_e(i));
    FF(2*i-1:2*i,1) = cell2mat(f_e(i));
    MM(2*i-1:2*i,2*i-1:2*i) = cell2mat(M_e(i));
end

%%build matrix 

L_g = A'*LL*A;
F_g = A'*FF;
M_g = A'*MM*A;


%%apply Direchlt BC in the code
L_gn = L_g(2:n,2:n);
F_gn = F_g(2:n);
F_gn(1) = F_gn(1) - L_g(2,1)*200;
F_gn(n-1) = F_gn(n-1) - L_g(n+1,n)*200;

M_gn = M_g(2:n,2:n);
M_gr1 = M_g(2:n,1);
M_grn = M_g(2:n,n+1);



 %%solution for poisson equation
                    T = 200+sparse(n+1,1);
                    T_g = L_gn\F_gn;
                    T(2:n,1) = T_g;
                    
                    
                    %%analytical
                    a = 100*(exp(L/2)-exp(-L/2))/L;
                    b = 200 + 100*exp(L/2) -a*L/2;
                    T_a = -100*exp(x) + a*x+b; 
                    
                    % plot(x,T,x,T_a,'--');
                    figure(3)
                    plot(x,T,x,T_a,'r--');
                    
                    legend('fem','analytical')


%%timeporal solution

N = M_g - (alpha*dt*L_g);
N_n = N(2:n,1:n+1);


for t = 0:dt:t_total
    T_temp = T;
    T0 = (N_n*T - 200*(M_gr1+M_grn));     %%I can't explain it in the code, just do it.
    T_n = M_gn\T0;
    T(2:n,1) = T_n;
    
    figure(3)
    plot(x,T,x,T_a,'r--');
    legend('fem','intial')
    ylim([200,220])
    drawnow
    t
    
    if norm(T - T_temp)/norm(T) < 10^-5         %%error has to be changed, smaller mesh, smaller tolerance
        break;
    end;
end
       
