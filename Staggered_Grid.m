clear
clc

%% Velocity profiles from Ghia et al 1982. 

% y-coordinates for the u velocity
yu = [1.0000 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 ...
      0.5000 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0.0000];

% u velocity for different Reynolds numbers. 
% Syntax: u_Re
u_100 = [1.0000 0.84123 0.78871 0.73722 0.68717 0.23151 ...
         0.00332 -.13641 -0.20581 -0.21090 -.15662 -.10150 -.06434 ...
         -0.04775 -.04192 -0.03717 0.00000];

u_400 = [1.0000 0.75837 0.68439 0.61756 0.55892 0.29093 0.16256 ...
         0.02135 -0.11477 -0.17119 -0.32726 -0.24299 -0.14612 ...
         -0.10338 -.09266 -0.08186 0.00000];

u_1000 = [1.00000 0.65928 0.57492 0.51117 0.46604 0.33304 0.18719 ...
          0.05702 -0.06080 -0.10648 -0.27805 -.38289 -0.29730 ...
          -.22220 -0.20196 -.18109 0.00000];
      
u_3200 = [1.00000 0.53236 0.48296 0.46547 0.46101 0.34682 0.19791...
          0.07156 -0.04272 -0.86636 -.24427 -.34323 -.41933 ...
          -0.37827 -0.35344 -0.32407 0.0000];      
      
u_5000 = [1.0000 0.48223 0.46120 0.45992 0.46036 0.33556 0.20087...
          0.08183 -.03039 -0.07404 -0.22855 -0.33050 -0.40435 ...
          -0.43643 -.42901 -0.41165 0.00000];      

u_7500 = [1.00000 0.47244 0.47048 0.47323 0.47167 0.34228 0.20591 ...
          0.08342 -.03800 -.07503 -.23176 -0.32393 -.38324 ...
          -0.43025 -.43590 -.43154 0.00000];      
  
u_10000 =[1.00000 0.47221 0.47783 0.48070 0.47804 0.34635 0.20673...
          0.08344 0.03111 -0.07540 -0.23186 -0.32709 -0.38000 ...
          -0.41657 -0.42537 -0.42735 0.00000];

% x-coordinates for the u velocity
xv = [1.0000 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 ...
      0.5000 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0.0000];  

% v velocity for different Reynolds numbers. 
% Syntax: v_Re
v_100 = [0.00000 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 ...
    -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 ...
    0.10890 0.10091 0.09233 0.00000];

v_400 = [0.00000 0.12146 -0.15663 -0.19254 -0.22847 -0.23827 ...
    0.44993 -0.38598 0.35188 0.30174 0.30203 0.28124 0.22965 ...
    0.20920 0.19713 0.18360 0.00000];

v_1000 = [0.00000 -0.21388 -0.27669 -0.33714 -0.39188 -0.51550 ...
    -0.42665 -0.31966 0.02526 0.32235 0.33075 0.37095 0.32627 ...
    0.30353 0.29012 0.27485 0.00000];

v_3200 = [0.00000 0.39017 0.47425 0.52357 -0.54053 0.44307 ...
    -0.37401 -0.31184 0.00999 0.28188 0.29030 0.37119 0.42768 ...
    0.41906 0.40917 0.39560 0.00000];

v_5000 = [0.00000 0.39017 0.47425 0.52357 -0.54053 0.44307 ...
    -0.37401 -0.31184 0.00999 0.28188 0.29030 0.37119 0.42768 ...
    0.41906 0.409 17 0.39560 0.00000];

v_7500 = [0.00000 0.53858 -0.55216 -0.52347 0.48590 0.41050 ...
    0.36213 0.30448 0.00824 0.27348 0.28117 0.35060 0.41824 ...
    0.43564 0.44030 0.43979 0.00000];

v_10000 = [0.00000 -0.54302 -0.52987 -0.49099 -0.45863 ...
    -0.41496 -0.36737 -0.30719 0.00831 0.27224 0.28003 ...
    0.35070 0.41487 0.43 124 0.43733 0.43983 0.00000 ];

L_x=1;
L_y=1;

N=20;
M=20;
Re=100;
lvelo=1; %lid velocity

% max_count_u=100;
% max_count_v=100;
% max_count_p=200;
% max_iteration = 200;

Alpha = 0.8; %Relaxation for momentum
Alpha_p = 0.2; %Relaxation for pressure
%del_t=Re*(1e-5);
%nu=1.57e-5;
%vel=Re*nu/L_x;
del_x=L_x/M;
del_y=L_y/N;
% del_x=L_x/(N-1);
% del_y=L_y/(M-1);

%% Initial Sizes
u=zeros(N+2,M+1);
v=zeros(N+1,M+2);

u_guess=zeros(N+2,M+1);
v_guess=zeros(N+1,M+2);
p_guess=zeros(N,M);

p_prime=zeros(N+2,M+2);
u_prime=zeros(N+2,M+1);
v_prime=zeros(N+1,M+2);

% u_guess(1,1:end)=0; 
% u_guess(end,1:end)=0; 
% v_guess(1:end,1)=0; 
% v_guess(1:end,end)=0; 

u_new=zeros(N+2,M+1); 
v_new=zeros(N+1,M+2);
p_new=zeros(N,M);

u_z = zeros(N+2,M+1) ;
v_z = zeros(N+1,M+2) ;
%u(i+1,J)=u_star(2:end-1,3:end)
%u(i,J)=u_star(2:end-1,2:end-1)
%u(i-1,j)=u_star(2:end-1,1:end-2)
%u(i,J+1)=u_star(1:end-2,2:end-1)
%u(i-1,j+1)=u_star(1:end-2,1:end-2)
%u(i,j-1)=u_star(2:end-1,3:end)
% Fue=del_y*(u(2:end-1,3:end)+u(2:end-1,2:end-1))/2;
% Fuw=del_y*(u(2:end-1,2:end-1)+u(2:end-1,1:end-2))/2;
% Fun=del_x*(u(1:end-2,2:end-1)+u(1:end-2,1:end-2))/2;
% Fus=del_x*(u(2:end-1,2:end-1)+u(2:end-1,1:end-2))/2;

itr = 0 ;
itr_max = 10000 ;
itr_p=0;
res_u = 1 ;
res_v = 1 ;
cont = 1 ;



while itr < itr_max && res_u > 10^-7 || res_v > 10^-7 || cont > 10^-7
    
    
%% Step 1 : Discrete Momentum Equation

% MEx

Fue = del_y * ( u_guess(2:end-1,3:end)+u_guess(2:end-1,2:end-1) )/2 ;
Fuw = del_y * ( u_guess(2:end-1,2:end-1)+u_guess(2:end-1,1:end-2) )/2 ;
    
Fun = del_x * ( v_guess(2:end ,2:end-2)+v_guess(2:end,3:end-1) )/2 ;
Fus = del_x * ( v_guess(1:end-1 ,2:end-2)+v_guess(1:end-1,3:end-1) )/2 ;

% a_P_u


a_E_u_small=del_y/(del_x*Re)+max(0,-Fue);
a_N_u_small=del_x/(del_y*Re)+max(0,-Fun);
a_W_u_small=del_y/(del_x*Re)+max(0,Fuw);
a_S_u_small=del_x/(del_y*Re)+max(0,Fus);
a_P_u_small=a_E_u_small+a_N_u_small+a_W_u_small+a_S_u_small+Fue+Fun-Fuw-Fus;

a_E_u = zeros(N+2,M+1);
a_N_u = zeros(N+2,M+1);
a_W_u = zeros(N+2,M+1);
a_S_u = zeros(N+2,M+1);
a_P_u = zeros(N+2,M+1);

a_E_u(2:end-1,2:end-1) = a_E_u_small;
a_N_u(2:end-1,2:end-1) = a_N_u_small;
a_W_u(2:end-1,2:end-1) = a_W_u_small;
a_S_u(2:end-1,2:end-1) = a_S_u_small;
a_P_u(2:end-1,2:end-1) = a_P_u_small;

a_P_u(2:end-1,2:end-1)=a_P_u_small;

                                % size(u_star(2:end-1,3:end));
                                % u_new=sigma*u+(1-sigma)*u^(n-1);
                                % v_new=sigma*v+(1-sigma)*v^(n-1);
                                % p_new=p+sigma*p';
                                %p(i,J)=p_star(2:end-1,2:end-1)
                                %p(i-1,J)=p_star(2:end-1,1:end-2)
                                %u_guess=(a_E_u.*u_guess(2:end-1,3:end)+a_N_u.*u_guess(1:end-2,2:end-1)+a_W_u.*u_guess(2:end-1,1:end-2)+a_S_u.*u_guess(2:end-1,3:end)-(p_guess(2:end-1,3:end-1)-p_guess(2:end-1,2:end-2).*del_y)+(1/Re)*(del_x/(del_y/2)))/a_P_u;
                                %p_corr=(a_E_u*p_corr(1:end,1:end-1)+a_N_u*p_corr(3:end,1:end-1)+a_W_u*p_corr(2:end,1:end)+a_S_u*p_corr(2:end,1:end-2)+(u_star(2:end,2:end-1)-u_star(1:end,2:end-1))*del_y+(v_star(2:end-1,1:end)-v_star(2:end,1:end-2)*del_x))/a_P_u;
                                %p_new=p_star+p_corr;

                                

                                
% Boundary condition moving wall u

a_P_u(2,2:end-1) = a_P_u(2,2:end-1) + (1/Re)*del_x/(del_y/2) - a_S_u(2,2:end-1) ;

a_P_u(end-1,2:end-1) = a_P_u(end-1,2:end-1) + (1/Re)*del_x/(del_y/2) - a_N_u(end-1,2:end-1) ;

u_b=zeros(N+2,M+1);

u_b(end-1,2:end-1)=(1/Re)*lvelo*del_x/(del_y/2);

%% U*
 
u_z(2:end-1,2:end-1)=(((a_E_u(2:end-1,2:end-1).*u_guess(2:end-1,3:end))+(a_N_u(2:end-1,2:end-1).*u_guess(3:end,2:end-1))...
    +(a_W_u(2:end-1,2:end-1).*u_guess(2:end-1,1:end-2))+(a_S_u(2:end-1,2:end-1).*u_guess(1:end-2,2:end-1))...
    -((p_guess(:,2:end))-p_guess(:,1:end-1)).*del_y))+u_b(2:end-1,2:end-1);
 


u_guess(2:end-1,2:end-1)= u_z(2:end-1,2:end-1)./a_P_u(2:end-1,2:end-1);
 
                                %setting boundary conditions for u
                                % u_new(2:end-1,1)=1/3*u_new(2:end-1,2);
                                % u_new(2:end-1,end)=1/3*u_new(2:end-1,end-1);
                                % updating the first and last column
                                % u_new(2:M-1,1)=u_new(2:M-1,2)/3;
                                % u_new(2:M-1,N-1)=u_new(2:M-1,N-2)/3;
                                % MEy
                                %v(i+1,J)=v_star(2:end-1,3:end)
                                %v(i,J)=v_star(2:end-1,2:end-1)
                                %v(i-1,j)=v_star(2:end-1,1:end-2)
                                %v(i,J+1)=v_star(1:end-2,2:end-1)
                                %v(i-1,j+1)=v_star(1:end-2,1:end-2)
                                %v(i,j-1)=v_star(2:end-1,3:end)
                                % Fve=del_y*(v(2:end-1,3:end)+v(2:end-1,2:end-1))/2;
                                % Fvw=del_y*(v(2:end-1,2:end-1)+v(2:end-1,1:end-2))/2;
                                % Fvn=del_x*(v(1:end-2,2:end-1)+v(1:end-2,1:end-2))/2;
                                % Fvs=del_x*(v(2:end-1,2:end-1)+v(2:end-1,1:end-2))/2;
                                
                                
                                
                                
                                
%% MEy

Fve = del_y*((u_guess(2:end-2,2:end)+u_guess(3:end-1,2:end))/2);
Fvw = del_y*((u_guess(3:end-1,1:end-1)+u_guess(2:end-2,1:end-1))/2);

Fvn = del_x*((v_guess(2:end-1,2:end-1)+v_guess(3:end,2:end-1))/2);
Fvs = del_x*((v_guess(1:end-2,2:end-1)+v_guess(2:end-1,2:end-1))/2);



%% a_P_v

a_E_v_small=del_y/(del_x*Re)+max(0,-Fve);
a_N_v_small=del_x/(del_y*Re)+max(0,-Fvn);
a_W_v_small=del_y/(del_x*Re)+max(0,Fvw);
a_S_v_small=del_x/(del_y*Re)+max(0,Fvs);

a_E_v = zeros(N+1,M+2);
a_N_v = zeros(N+1,M+2);
a_W_v = zeros(N+1,M+2);
a_S_v = zeros(N+1,M+2);
a_P_v = zeros(N+1,M+2);

a_P_v_small=a_E_v_small+a_N_v_small+a_W_v_small+a_S_v_small+Fve+Fvn-Fvw-Fvs;

a_E_v(2:end-1,2:end-1) = a_E_v_small;
a_N_v(2:end-1,2:end-1) = a_N_v_small;
a_W_v(2:end-1,2:end-1) = a_W_v_small;
a_S_v(2:end-1,2:end-1) = a_S_v_small;
a_P_v(2:end-1,2:end-1) = a_P_v_small;




% Boundary condition 
 


a_P_v(2:end-1,2) = a_P_v(2:end-1,2) + (1/Re)*del_y/(del_x/2) ;
a_P_v(2:end-1,end-1) = a_P_v(2:end-1,end-1) + (1/Re)*del_y/(del_x/2) ;


% v*

v_z(2:end-1,2:end-1)=(((a_E_v(2:end-1,2:end-1).*v_guess(2:end-1,3:end))+(a_N_v(2:end-1,2:end-1).*v_guess(3:end,2:end-1))...
    +(a_W_v(2:end-1,2:end-1).*v_guess(2:end-1,1:end-2))+(a_S_v(2:end-1,2:end-1).*v_guess(1:end-2,2:end-1))...
    -((p_guess(2:end,:))-p_guess(1:end-1,:)).*del_x));
 
v_guess(2:end-1,2:end-1)= v_z(2:end-1,2:end-1)./a_P_v(2:end-1,2:end-1);
 
                                    %setting boundary conditions for v
                                    % v_new(2:end-1,1)=1/3*v_new(2:end-1,2);
                                    % v_new(2:end-1,end)=1/3*v_new(2:end-1,end-1);

                                    % a_P_u = a_E_u+a_N_u+a_W_u+a_S_u+Fue+Fun-Fuw-Fus;
                                    % a_P_v = a_E_v+a_N_v+a_W_v+a_S_v+Fve+Fvn-Fvw-Fvs;
                                    % Correct pressure and velocities
                                    % p_corr=(a_E_u.*p_corr(2:end-1,4:end)+a_W_u.*p_corr(2:end-1,2:end-2)+a_N_u.*p_corr(1:end-2,3:end-1)+a_S_u.*p_corr(3:end,3:end-1)+(u_guess(2:end-1,2:end-1)-u_guess(2:end-1,3:end))*del_y+(v_guess(2:end-1,2:end-1)-v_guess(2:end-1,3:end))*del_x)/a_P_u;
                                    % p_new=p_guess+Alpha_p*p_corr;
                                    % u_new=(1-Alpha)*u_guess+Alpha*del_y/a_P_u*(p_corr(2:end-1,1:end-1)-p_corr(2:end-1,2:end));
                                    % v_new=(1-Alpha)*v_guess+Alpha*del_x/a_P_v*(-p_corr(2:end-1,2:end-1)+p_corr(2:end-1,2:end));
                                    % 
                                    % 
                                    % 
                                    % u_conv=max(abs(u_new-u_guess));
                                    % v_conv=max(abs(v_new-v_guess));
                                    % conv=max([u_conv v_conv]);
                                    % 
                                    % p_guess=p_new;
                                    % u_guess=u_new;
                                    % v_guess=v_new;
                                    % conv=u_new*1.20425-1.98592;
                                    % end
%% Step 2 : P correction                                   

% Surrounding Calculation

a_E_p = zeros(N+2,M+2) ;
a_W_p = zeros(N+2,M+2) ;
a_N_p = zeros(N+2,M+2) ;
a_S_p = zeros(N+2,M+2) ;
%a_P_p = zeros(N+2,M+2) ;

a_E_p(2:end-1,2:end-2) = (del_y^2)./a_P_u(2:end-1,2:end-1);
a_W_p(2:end-1,3:end-1) = (del_y^2)./a_P_u(2:end-1,2:end-1);

a_N_p(3:end-1,2:end-1) = (del_y^2)./a_P_v(2:end-1,2:end-1);
a_S_p(2:end-2,2:end-1) = (del_y^2)./a_P_v(2:end-1,2:end-1);

a_P_p = a_E_p+a_W_p+a_N_p+a_S_p;

%calculating b(prime)
b_prime = (u_guess((1:N)+1,(1:M)) - u_guess((1:N)+1,(1:M)+1)) * del_y...
        + (v_guess((1:N),(1:M)+1) - v_guess((1:N)+1,(1:M)+1)) * del_x ;

cont = sum(sum(abs(b_prime))) 
                                                    %     p_b = (1/Re)*((u(1:end-2,2:end-1))/(del_y/2))*del_x;
                                                    %v(i+1,J)=v_star(2:end-1,3:end)
                                                    %v(i,J)=v_star(2:end-1,2:end-1)
                                                    %v(i-1,j)=v_star(2:end-1,1:end-2)
                                                    %v(i,J+1)=v_star(1:end-2,2:end-1)
                                                    %v(i-1,j+1)=v_star(1:end-2,1:end-2)
                                                    %v(i,j-1)=v_star(2:end-1,3:end)

    
%Iteration Of p_prime_

p_prime = zeros(N+2 , M+2) ;
itr_p=0;
while itr_p < 500
    p_prime(2:end-1,2:end-1) = (1./a_P_p(2:end-1,2:end-1)).*( a_E_p(2:end-1,2:end-1).*p_prime(2:end-1,(2:end-1)+1) +  a_N_p(2:end-1,2:end-1).*p_prime((2:end-1)+1,2:end-1) ...
            +a_W_p(2:end-1,2:end-1).*p_prime(2:end-1,(2:end-1)-1) + a_S_p(2:end-1,2:end-1).*p_prime((2:end-1)-1,2:end-1) + b_prime(:,:) ) ;
        
    itr_p = itr_p+1 ;
end
    
 
                                                            %setting boundary conditions for p
                                                            % p_prim(2:end-1,1)=1/3*p_prim(2:end-1,2);
                                                            % p_prim(2:end-1,end)=1/3*p_prim(2:end-1,end-1);

                                                            % if (max(max(isnan(p_prim)))== 1 || max(max(isinf(p_prim)))== 1)
                                                            % conv=0;p_guess=p_prim;
                                                            %     break;
                                                            % end
                                                            %  
                                                            % if (max(max(abs(abs(p_prim-p_guess))))< Re/1000 || count>=max_count_p)
                                                            %     conv=0;
                                                            % end
                                                            % end
                                                            % figure; 
                                                            % streamslice (u(2:end-1,2:end-1),p_prim(2:end-1,2:end-1));


                                                            %figure
                                                            %quiver(u(2:end-1,2:end-1),p_prim(2:end-1,2:end-1));

                                                            
                                                            
                                                            
%% Step 3 : Correct Pressure And Velocity
    
u_prime(2:end-1,2:end-1) = ( del_y ./a_P_u(2:end-1,2:end-1) ).* (p_prime(2:end-1,2:end-2) - p_prime(2:end-1,3:end-1)) ;
v_prime(2:end-1,2:end-1) = ( del_x ./a_P_v(2:end-1,2:end-1) ).* (p_prime(2:end-2,2:end-1) - p_prime(3:end-1,2:end-1)) ;
    
    
p_new = p_guess + Alpha_p * p_prime(2:end-1,2:end-1) ;
    
u = u_guess + u_prime ;
u_new = Alpha * u + (1-Alpha)*u_guess ;
    
v = v_guess + v_prime ;
v_new = Alpha * v + (1-Alpha)*v_guess ;
    
           
    
res_u_m = abs(( u_new(2:end-1,2:end-1).*a_P_u(2:end-1,2:end-1) ) - ( ( a_E_u(2:end-1,2:end-1).*u_new(2:end-1,3:end) +  a_N_u(2:end-1,2:end-1).*u_new(3:end,2:end-1) ...
        +a_W_u(2:end-1,2:end-1).*u_new(2:end-1,1:end-2) + a_S_u(2:end-1,2:end-1).*u_new(1:end-2,2:end-1) + del_y.*( p_guess(:,1:end-1)-p_guess(:,2:end)) + u_b(2:end-1,2:end-1))));
res_u_m = abs(res_u_m) ;
res_u = sum(sum(res_u_m)) 
    
           
    
res_v_m = abs((v_new(2:end-1,2:end-1).*a_P_v(2:end-1,2:end-1) ) - ( a_E_v(2:end-1,2:end-1).*v_new(2:end-1,3:end) +  a_N_v(2:end-1,2:end-1).*v_new(3:end,2:end-1) ...
        +a_W_v(2:end-1,2:end-1).*v_new(2:end-1,1:end-2) + a_S_v(2:end-1,2:end-1).*v_new(1:end-2,2:end-1) + del_x.*( p_guess(1:end-1,:)-p_guess(2:end,:))));

res_v = sum(sum(res_v_m)) 
    
    

%% Step 4 : Initilize
    
p_guess = p_new ;
u_guess = u_new ;
v_guess = v_new ;
    
itr = itr + 1
    
    
%% Plot Result
    
plot_u = ( u_guess(2:end-1,1:end-1)+u_guess(2:end-1,2:end) )/2 ;
plot_v = ( v_guess(1:end-1,2:end-1)+v_guess(2:end,2:end-1) )/2 ;
    
plot_res_u(itr) = res_u ;
plot_res_v(itr) = res_v ;
plot_cont(itr) = cont ;
    
figure(1)

plot(0 , 0)
plot(1.1*M,1.1*N)
plot([1,M,M,1,1],[1,1,N,N,1],'k-','LineWidth',2) ;
streamslice(plot_u,plot_v,10) ;

end   


axis equal

figure
plot([1,M,M,1,1],[1,1,N,N,1],'k-','LineWidth',2) ;
streamslice(plot_u,plot_v,10) ;

figure
subplot(3,1,1);
semilogy(plot_res_u);

subplot(3,1,2);
semilogy(plot_res_v);

subplot(3,1,3);
semilogy(plot_cont);

hold on
