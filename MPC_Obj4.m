%Project, Name: Anshumali Jaiswal, Roll nubmer: 213020033
function obj  = MPC_Obj4(Uf)

global p q Wy Wdu nst nop nip
global xk_hat rkf b_hat C_mat ukminus_1
global sys

zk_hat = zeros(nst,p);
zk_hat(:,1) = xk_hat;
yk_hat = zeros(nop,p);
yk_hat(:,1) = C_mat*zk_hat(:,1);

uf = zeros(nip,p);
n2 = 0;
for k = 1:p
    if k <= q
        n1 = n2 + 1;
        n2 = k*nip;
        uf(:,k) = Uf(n1:n2,1)';
    else
        uf(:,k) = Uf(n1:n2,1)';
    end
end

deluk = uf(:,1) - ukminus_1;
err = rkf- yk_hat(:,1);
obj = err'*Wy*err + deluk'*Wdu*deluk;  %Objective function
for k = 2:p
    zk_hat(:,k) = sys.phy*zk_hat(:,k-1) + sys.gama_u*(uf(:,k-1) + b_hat);
    yk_hat(:,k) = C_mat*zk_hat(:,k);
    err = rkf - yk_hat(:,k);
    deluk = uf(:,k) - uf(:,k-1);
    obj = obj + err'*Wy*err + deluk'*Wdu*deluk;
end

