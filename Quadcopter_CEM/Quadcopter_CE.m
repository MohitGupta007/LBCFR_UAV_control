%% ============================================================ %%
%%  Quadcopter MPPI (Model Predictive Path Integral Control)
%% ============================================================ %%

function [u_seq,sigma2new] = Quadcopter_CE(state,model,params,previous_u_seq)

% ======================================================================= %
%   General Parameter Setting
% ======================================================================= %

K  = params.K;
Hp = params.Hp;
control_dt = params.control_dt;
plant_dt = params.plant_dt;
num_elites = 5;
nu = 1000;
R  = 0.01;
num_iterations = params.num_iterations;

state_dim = size(state,1);
u_dim = 4;

u_init = zeros(u_dim,1);
sigma = params.sigma;

switch nargin
    case 3
        u_seq = 620.6108*ones(u_dim,Hp);
    case 4
        u_seq(:,1:Hp-1) = previous_u_seq(:,2:Hp);
        u_seq(:,Hp) = previous_u_seq(:,Hp);
        sigma(:,1:Hp-1) = params.sigma(:,2:Hp);
end
% sigma = params.sigma;
S_seq = zeros(K,Hp);

x_seq = zeros(state_dim,K,Hp);
state0 = repmat(state,1,K);
x_seq(:,:,1) = state0;
x_plant_seq = zeros(state_dim,K,(control_dt/plant_dt)+1);
du_seq = zeros([u_dim,K,Hp]);
udu = zeros(u_dim,K,Hp);

for it = 1:num_iterations
    
for i=1:Hp
    dT_seq = normrnd(0,real(sqrt(sigma(1,i))),[1,K]);
    dTau_phi_seq = normrnd(0,real(sqrt(sigma(2,i))),[1,K]);
    dTau_theta_seq = normrnd(0,real(sqrt(sigma(3,i))),[1,K]);
    dTau_psi_seq = normrnd(0,real(sqrt(sigma(4,i))),[1,K]);

%     du_seq(1,:,i) = dT_seq - dTau_theta_seq - dTau_psi_seq;
%     du_seq(2,:,i) = dT_seq - dTau_phi_seq + dTau_psi_seq;
%     du_seq(3,:,i) = dT_seq + dTau_theta_seq - dTau_psi_seq;
%     du_seq(4,:,i) = dT_seq + dTau_phi_seq + dTau_psi_seq;
    
    du_seq(1,:,i) = dT_seq;
    du_seq(2,:,i) = dTau_phi_seq;
    du_seq(3,:,i) = dTau_theta_seq;
    du_seq(4,:,i) = dTau_psi_seq;
    
    udu(:,:,i) = u_seq(:,i)+du_seq(:,:,i);

    x_plant_seq(:,:,1) = x_seq(:,:,i);
    for j=1:(control_dt/plant_dt)
        x_plant_seq(:,:,j+1) = get_new_states(x_plant_seq(:,:,j),u_seq(:,i)+du_seq(:,:,i),plant_dt,model);
    end
    x_seq(:,:,i+1) = x_plant_seq(:,:,end);
end

for k=1:K
    S_temp_seq = zeros(1,Hp);
    S_temp_seq(1,Hp) = Quadcopter_Costfunc(x_seq(:,k,Hp),u_seq(:,Hp),du_seq(:,k,Hp),nu,R,params.xd(Hp,:));

    for i = 1:Hp-1
        S_temp_seq(1,Hp-i) = S_temp_seq(1,Hp+1-i) + Quadcopter_Costfunc(x_seq(:,k,Hp-i),u_seq(:,Hp-i),du_seq(:,k,Hp-i),nu,R,params.xd(Hp-i,:));
    end

    S_seq(k,:) = S_temp_seq;    
end
J = S_seq(:,1);
[~,Ie] = sort(J);
Ie = Ie(1:num_elites);
Je= J(Ie);  
%find the sample with the minimum cost
[Jmin, Imin] = min(J);    
%select the parameters with minimum cost
Umin = udu(:,Imin,:); 
%select the elite parameters (FILL IN)
udue = udu(:,Ie,:);

%use maximum likelihood esitmation to compute the new parameter
%distribution (FILL IN)
for j=1:Hp
    Unew(:,j) = mean(udue(:,:,j),2);
    sigma2new(:,j) = sum((udue(:,:,j)-Unew(:,j)).^2,2)./(num_elites) + 0.001;
end
% u_new = Unew(:,1)'
% rollout the mean for the new parameter distribution and evaluate the
% cost
 for i = 1:u_dim
    u_seq(i,:) = real(normrnd(0,real(sqrt(sigma2new(i,:))),1,Hp) + Unew(i,:));

    u_max_index = (u_seq(i,:) > 900);
    u_seq(i,u_max_index) = 900;

    %u_seq(i,:) = u_seq(i,:) + du_seq(S_min_index,:,i);
 end
sigma = sigma2new;
end

end
