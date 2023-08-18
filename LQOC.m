%Project, Name: Anshumali Jaiswal, Roll nubmer: 213020033 

close all
clear all
clc

global sys
load System2_Parameters

load ARXLQOCParameters

delta_t  = 0.1; %sampling time in mins
tf       = 25 ; %simulation time
N        = 250; %No of samples
randn('seed',0);

nst = 3;  %no of states
nip = 2;  %no of manipulated inputs
nd  = 1;  %no of disturbances
nop = 2;  %no of outputs

Xk_state = zeros(nst,N);   %state profiles
Uk_man   = zeros(nip,N);   %input profiles
Dk_dis   = zeros(nd,N);    %disturbance profile
Yk_out   = zeros(nop,N);   %output profile
Rk       = zeros(nop,N);
rk       = zeros(nop,N);
ek       = zeros(nop,N);

R_mat = diag(sys.meas_sigma.^2);

%Initialisation
Xk_state(:,1) = sys.Xs + [0.03 5 0.01]';
Uk_man(:,1)   = sys.Us;
Dk_dis(:,1)   = sys.Ds;
Yk_out(:,1)   = sys.Ys;
z_hat(:,1)    = zeros(6,1);

%Upper and Lower limits
Uk_H = [2 430]';
Uk_L = [0 350]';
uk_H = Uk_H - sys.Us;
uk_L = Uk_L - sys.Us;


%dynamic Simulation

Tj      = zeros(1,N);  %to record time
uk(:,1) = zeros(nip,1);
vk      = zeros(nop,1);

v  = [1 100];
Wy = diag(v);

b   = [1 10];
eta = 0.4;
Wu  = eta*diag(b);

Ku             = Cmimo*inv(eye(6) - phyMimo)*gammaMimo;
Ke             = Cmimo*inv(eye(6) - phyMimo)*L + eye(2);
[Ginf,Sinf,Ev] = dlqr(phyMimo,gammaMimo,(Cmimo'*Wy*Cmimo),Wu);
alpha          = 0.9;

for k = 1:N-1
    Tj(k) = (k-1)*delta_t; %update time
       
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

    d  = (mvnrnd(0,sys.dk_sigma^2,1))';
    Dk_dis(k) = sys.Ds + dstep + d ;

    Rk(:,k) = sys.Ys + rk(:,k);

    y(:,k) = Yk_out(:,k) - sys.Ys;

    ek(:,k) = y(:,k) - Cmimo*z_hat(:,k);
    phi_e = alpha*eye(2);

    if k == 1
       e_f(:,k) = (eye(2) - phi_e)*ek(:,k);
    else
       e_f(:,k) = phi_e*e_f(:,k-1) + (eye(2) - phi_e)*ek(:,k);
    end

    us(:,k) = inv(Ku)*(rk(:,k) - Ke*e_f(:,k));
    zs(:,k) = inv(eye(6) - phyMimo)*(gammaMimo*us(:,k) + L*e_f(:,k));

    uk(:,k) = us(:,k) - Ginf*(z_hat(:,k) - zs(:,k));

    if (uk(1,k)<= uk_L(1))
       uk(1,k) = uk_L(1);
    elseif(uk(1,k)>=uk_H(1))
       uk(1,k) = uk_H(1);
    end
    if (uk(2,k)<= uk_L(2))
       uk(2,k) = uk_L(2);
    elseif(uk(2,k)>=uk_H(2))
       uk(2,k) = uk_H(2);
    end
    
    Uk_man(:,k) = uk(:,k) + sys.Us;
    
    sys.Uk = Uk_man(:,k);  
    sys.Dk = Dk_dis(:,k);

    [t, Xt] = ode45('System2_Dynamics',[0 delta_t],Xk_state(:,k));
   
    Xk_state(:,k+1) = Xt(end,:);
   
    %noise in measurement

    vk = (mvnrnd(zeros(nop,1),R_mat,1))';
    Yk_out(:,k+1) = sys.C_mat*Xk_state(:,k+1) + vk;
     

    z_hat(:,k+1) = phyMimo*z_hat(:,k) + gammaMimo*uk(:,k) + L*ek(:,k);
    
    
end

%input at the final instant
Tj(N) = N*delta_t;
Uk_man(:,N) = Uk_man(:,N-1);
dk = sys.dk_sigma*randn;
Dk_dis(N) = Dk_dis(N-1);
uk(:,N) = uk(:,N-1);
us(:,N) = us(:,N-1);
Yk_out(:,N) = Yk_out(:,N-1);
Rk(:,N)     = Rk(:,N-1);


SSE = [0,0]';
for i = 1:250 
    SSE = SSE + (Yk_out(:,i) -Rk(:,i)).^2;
end
SSMV = [0,0]';
for j =1:250
    SSMV = SSMV + (uk(:,j)).^2;
end

%Display of results
for i =1:2
    figure(i)
    plot(Tj,Yk_out(i,:),'b-','LineWidth',1.5),grid,xlabel('Time(min)'),title('Output Trajectories')
    hold on
    plot(Tj,Rk(i,:),'r-','LineWidth',1.5),,legend('YLQOC','Setpoint')

    if i == 1
        ylabel('Y_1 and Rk_1')
    else
        ylabel('Y_2 and Rk_2')
    end
    hold off
end

figure(3),stairs(Tj,Dk_dis(1,:),'g-','LineWidth',1.5),grid on, ylabel('Dk'), xlabel('Time(min)'), title('Disturbances')
figure(4)
subplot(2,1,1),plot(Tj,Uk_man(1,:),'r-','LineWidth',1.5),grid on, ylabel('U_1k'), xlabel('Time(min)'), title('Input Trejactories')
subplot(2,1,2),plot(Tj,Uk_man(2,:),'m-','LineWidth',1.5),grid on, ylabel('U_2k'), xlabel('Time(min)'), title('Input Trejactories')