%% ============================================================ %%
%%  Quadcopter MPPI (Model Predictive Path Integral Control)
%% ============================================================ %%

function u_seq = Quadcopter_MPPI(state,model,params,previous_u_seq)

% ======================================================================= %
%   General Parameter Setting
% ======================================================================= %

K  = params.K;
Hp = params.Hp;
control_dt = params.control_dt;
plant_dt = params.plant_dt;
nu = 1000;
R  = 0.01;
lambda = 1000;

state_dim = size(state,1);
u_dim = 4;

% ======================================================================= %
%   Generate random control variations >>> (N*K)
% ======================================================================= %

%disturbance_parameter = 1;
%epsilon_sigma = 0.5;

if (state(1)-5)^2+(state(2)-5)^2+(state(3)-15)^2 <= 2
    dT_seq = normrnd(0,120,[K,Hp]);
    dTau_phi_seq = normrnd(0,10,[K,Hp]);
    dTau_theta_seq = normrnd(0,10,[K,Hp]);
    dTau_psi_seq = normrnd(0,10,[K,Hp]);
else
    dT_seq = normrnd(0,120,[K,Hp]);
    dTau_phi_seq = normrnd(0,10,[K,Hp]);
    dTau_theta_seq = normrnd(0,10,[K,Hp]);
    dTau_psi_seq = normrnd(0,10,[K,Hp]);
end
    
du_seq = zeros([u_dim,K,Hp]);
du_seq(1,:,:) = dT_seq - dTau_theta_seq - dTau_psi_seq;
du_seq(2,:,:) = dT_seq - dTau_phi_seq + dTau_psi_seq;
du_seq(3,:,:) = dT_seq + dTau_theta_seq - dTau_psi_seq;
du_seq(4,:,:) = dT_seq + dTau_phi_seq + dTau_psi_seq;
%du_seq = disturbance_parameter*1/sqrt(dt)*normrnd(0,epsilon_sigma,[u_dim,K,Hp]);

% ======================================================================= %
%   The cost of a trajectory
% ======================================================================= %

u_init = zeros(u_dim,1);

switch nargin
    case 3
        u_seq = 620.6108*ones(u_dim,Hp);
    case 4
        u_seq(:,1:Hp-1) = previous_u_seq(:,2:Hp);
        u_seq(:,Hp) = previous_u_seq(1,Hp);
end

S_seq = zeros(K,Hp);

x_seq = zeros(state_dim,K,Hp);
state0 = repmat(state,1,K);
x_seq(:,:,1) = state0;
x_plant_seq = zeros(state_dim,K,(control_dt/plant_dt)+1);

for i=1:Hp-1
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

du_seq = permute(du_seq,[2 3 1]);
[S_min_value, S_min_index] = min(S_seq(:,1));

S_min_seq = min(S_seq);
S_mean_seq = mean(S_seq);
S_std_seq = std(S_seq);
%S_seq = S_seq-S_min_seq+1;
for i = 1:u_dim
    u_seq(i,:) = u_seq(i,:) + (sum((exp(-(1/lambda)*S_seq)).*du_seq(:,:,i))./sum(exp(-(1/lambda)*S_seq)));
    
    u_max_index = (u_seq(i,:) > 900);
    u_seq(i,u_max_index) = 900;
    
    %u_seq(i,:) = u_seq(i,:) + du_seq(S_min_index,:,i);
end

end
