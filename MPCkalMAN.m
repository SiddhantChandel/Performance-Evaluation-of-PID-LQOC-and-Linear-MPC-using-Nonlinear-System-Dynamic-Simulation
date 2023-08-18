%Project, Name: Anshumali Jaiswal, Roll nubmer: 213020033

clear all
clc
close all

global p q Wy Wdu nst nop nip
global xk_hat rkf b_hat C_mat rk ukminus_1
global sys

load System2_Parameters.mat
load System2_Continuous_LinMod.mat
nst = 3; %no of states
nip = 2; %no of manipulated inputs
nd  = 1; %no of disturbances
nop = 2; %no of outputs

delta_t  = 0.1; %sampling time in mins
N        = 250;
Tj(N)        = 25; %in mins
randn('seed',0);
sys.phy    = expm(A_mat*delta_t);
sys.gama_u = (sys.phy - eye(nst))*inv(A_mat)*B_mat;
sys.gama_d = (sys.phy - eye(nst))*inv(A_mat)*H_mat;
sys.C_mat  = C_mat;

% Pole placement controller design 

ctrb_mat  =  ctrb(sys.phy,sys.gama_u);
n_ctrb    =  rank(ctrb_mat) ;
fprintf('\n Rank of controllability matrix: %d \n',n_ctrb);
con_poles =  [0.82 0.9 0.92]';
G_mat     =  place(sys.phy,sys.gama_u,con_poles);
eig(sys.phy-sys.gama_u*G_mat);

% state argumention for input bias formulation

gama_beta =  sys.gama_u;
phi_aug   =  [sys.phy  gama_beta ; zeros(nip,nst)  eye(nip)];
Cmat_aug  =  [sys.C_mat zeros(nop,nip)];
gama_aug  =  [sys.gama_u  ;  zeros(nip,nip)];
Ku        =  sys.C_mat*(inv(eye(3)-sys.phy))*sys.gama_u;
obsv_mat  =  obsv(phi_aug,Cmat_aug);                      % construct of observability matrix 
n_obsv    =  rank(obsv_mat);                              % rank of obsevability matrix 
fprintf('\n Rank of argumented observability matrix : %d \n ',n_obsv);

obs_poles =  [0.3 0.4 0.45 0.5 0.6]';
L_mat_aug =  place(phi_aug',Cmat_aug',obs_poles);
L_mat_aug =  L_mat_aug';
eig(phi_aug-L_mat_aug*Cmat_aug);


Xk_state = zeros(nst,N); %state profiles
Uk_man   = zeros(nip,N); %input profiles
Dk_dis   = zeros(nd,N);  %disturbance profile
Yk_out   = zeros(nop,N);  %output profile

uk = zeros(nip,N);
yk = zeros(nop,N);

Rk = zeros(nop,N);
rk = zeros(nop,N);
ek = zeros(nop,N);

Uk_H = [2 430]';
Uk_L = [0 350]';
uk_H = Uk_H - sys.Us;
uk_L = Uk_L - sys.Us;


% Initialisation

Xk_state(:,1) = sys.Xs + [0.03 5 0.01]';
Uk_man(:,1)   = sys.Us;
Dk_dis(:,1)   = sys.Ds;
Yk_out(:,1)   = sys.Ys;
x_aug(:,1)    = [0 0 0 0 0]';
d_aug(:,1)    = [0 0 0]';
Qd            = 0.015^2;
Q_beta        = (2/3)*[0.001^2  0 ; 0  0.05^2];

%dynamic Simulation
Tj = zeros(1,N);  %to record time

gamad_aug  = [sys.gama_d  zeros(3,2); zeros(2,1) eye(2)];
D_aug      = zeros(2,5);
Q_aug      = [Qd zeros(1,2);zeros(2,1) Q_beta];

N_aug = zeros(3,2);
R     = diag(sys.meas_sigma.^2);

dmod = ss(phi_aug,[gama_aug,gamad_aug],Cmat_aug,D_aug,0.1);

[KEST,L_aug, Pinf_aug] = kalman(dmod,Q_aug,R,N_aug);

% MPC Tuning Parameters
a  = [1 100];
Wy =  diag(a);
p  = 40;    %prediction horizon
q  = 5;     %control horizon
ukminus_1 = zeros(nip,1);
b   = [5 5];
Wdu = 0.4*diag(b);  %input move weighting matrix

% Inputconstraints
Uf_min   = zeros(q*nip,1);   %Future Input vector
U_min    = [];               %Future input lower bound vector 
U_max    = [];               %Future input upper bound vector 
    

for i = 1:q
    U_min    = [U_min;uk_L];

    U_max    = [U_max ; uk_H];

end


% Optimisation parameters
oldopts   = optimset;
tolerance = 10^(-8);
options   = optimset(oldopts,'MaxIter',1e6,'Display','iter','TolFun',tolerance,'LargeScale','off','TolX',tolerance);

rkf     = zeros(nop,N);
Rk(:,1) = rkf(:,1) + sys.Ys;

for k = 1:N-1
   Tj(k) = (k-1)*delta_t; %update time
     
   %Set point and disturbances
   if(k>=60 && k<=120)
           rk(:,k) = [10 -0.025]';
   else
           rk(:,k) = [0 0]';
   end
    
   if (k>=180)
           dstep = 0.5;
   else
           dstep = 0;
   end

   b_hat = x_aug(4:5,k);
   xk_hat = x_aug(1:3,k);
   
   if (k > 1)
       ukminus_1 = uk(:,k-1);
   end
   
   rkf = rk(:,k);

   n1 = nip + 1;
   n2 = (q-1)*nip + 1;
   Ufk_0 = [Uf_min(n1:end) ; Uf_min(n2:end)];
   Uf_min = fmincon('MPC_Obj4',Ufk_0,[],[],[],[],U_min,U_max,[],options);
   uk(:,k) = Uf_min(1:nip);
   Uk_man(:,k) = sys.Us + uk(:,k);
   
   d = (mvnrnd(0,sys.dk_sigma^2,1))';
   Rk(:,k) = sys.Ys + rk(:,k);
   Dk_dis(k) = sys.Ds + dstep + d ;

   
   
   yk(:,k) = Yk_out(:,k) - sys.Ys;
   ek(:,k) = yk(:,k) - Cmat_aug*x_aug(:,k); %innovation

   us(:,k) = inv(Ku)*rk(:,k) - x_aug(4:5,k);
   xs(:,k) = inv(eye(3) - sys.phy)*sys.gama_u*inv(Ku)*rk(:,k);

  
   
   sys.Uk = Uk_man(:,k);  
   sys.Dk = Dk_dis(:,k);
   
   [t, Xt] = ode45('System2_Dynamics',[0 delta_t],Xk_state(:,k));
   
   %Simulate plant dynamics
   Xk_state(:,k+1) = Xt(end,:);
   vk = (mvnrnd(zeros(nop,1),R,1))';
   Yk_out(:,k+1) = sys.C_mat*Xk_state(:,k+1) + vk; 
   
   %Kalman Predictor
   x_aug(:,k+1) = phi_aug*x_aug(:,k) + gama_aug*uk(:,k) + L_aug*ek(:,k);

 

end

% us(:,N) = us(:,N-1);
Dk_dis(:,N) = Dk_dis(:,N-1);
Yk_out(:,N) = Yk_out(:,N-1);
Uk_man(:,N) = Uk_man(:,N-1);
Rk(:,N)     = Rk(:,N-1);
Tj(N)       = N*delta_t;

%Performance indices

SSE = [0,0]';
for i = 1:250 
    SSE = SSE + (Yk_out(:,i) - Rk(:,i)).^2;
end
SSMV = [0,0]';
for j =1:250
    SSMV = SSMV + (uk(:,j)).^2;
end

% file='q=40';
% save(file,'Yk_out')

%Display of results
for i =1:2
    figure(i)
    plot(Tj,Yk_out(i,:),'b-','LineWidth',1.5),grid,xlabel('Time(min)'),title('Output Trajectories')
    hold on
    plot(Tj,Rk(i,:),'r-','LineWidth',1.5),legend('YMPC','Setpoint')

    if i == 1
        ylabel('Y_1 and Rk_1')
    else
        ylabel('Y_2 and Rk_2')
    end
    hold off
end

figure(3),stairs(Tj,Dk_dis,'g-','LineWidth',1.5),grid on, ylabel('Dk'), xlabel('Time(min)'), title('Disturbances')
figure(4)
subplot(2,1,1),plot(Tj,Uk_man(1,:),'r-','LineWidth',1.5),grid on, ylabel('U_1k'), xlabel('Time(min)'), title('Input Trejactories')
subplot(2,1,2),plot(Tj,Uk_man(2,:),'m-','LineWidth',1.5),grid on, ylabel('U_2k'), xlabel('Time(min)'), title('Input Trejactories')