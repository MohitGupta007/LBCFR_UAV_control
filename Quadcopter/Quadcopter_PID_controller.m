%% ============================================================ %%
%%  Simple PID Controller
%% ============================================================ %%

function [control, T_and_Tau] = Quadcopter_PID_controller(state,state_diff,pre_control)

g = 9.81;
m  = 0.468;
L  = 0.225;
kF = 2.980*10^(-6);
kM = 1.140*10^(-7);
Ixx = 4.856*10^(-3);
Iyy = 4.856*10^(-3);
Izz = 8.801*10^(-3);

x       = state(1);
y       = state(2);
z       = state(3);
vx      = state(4);
vy      = state(5);
vz      = state(6);
phi     = state(7);
theta   = state(8);
psi     = state(9);
p       = state(10);
q       = state(11);
r       = state(12);

x_dif       = state_diff(1);
y_dif       = state_diff(2);
z_dif       = state_diff(3);
vx_dif      = state_diff(4);
vy_dif      = state_diff(5);
vz_dif      = state_diff(6);
phi_dif     = state_diff(7);
theta_dif   = state_diff(8);
psi_dif     = state_diff(9);
p_dif       = state_diff(10);
q_dif       = state_diff(11);
r_dif       = state_diff(12);

K_z_D      = 2.5;
K_phi_D    = 1.75;
K_theta_D  = 1.75;
K_psi_D    = 1.75;
K_z_P      = 1.5;
K_phi_P    = 6;
K_theta_P  = 6;
K_psi_P    = 6;

z_d = 20;
phi_d = 0;
theta_d = 0;
psi_d = 0;

T         = (g + K_z_D*(0-z_dif) + K_z_P*(z_d-z))*m/cos(phi)/cos(theta);
Tau_phi   = (K_phi_D*(0-phi_dif) + K_phi_P*(phi_d-phi))*Ixx;
Tau_theta = (K_theta_D*(0-theta_dif) + K_theta_P*(theta_d-theta))*Iyy;
Tau_psi   = (K_psi_D*(0-psi_dif) + K_psi_P*(psi_d-psi))*Izz;

control(1) = T/(4*kF) - Tau_theta/(2*kF*L) - Tau_psi/(4*kM);
control(2) = T/(4*kF) - Tau_phi/(2*kF*L) + Tau_psi/(4*kM);
control(3) = T/(4*kF) + Tau_theta/(2*kF*L) - Tau_psi/(4*kM);
control(4) = T/(4*kF) + Tau_phi/(2*kF*L) + Tau_psi/(4*kM);

T_and_Tau = sqrt([abs(T/(4*kF));abs(Tau_phi/(2*kF*L));abs(Tau_theta/(2*kF*L));abs(Tau_psi/(4*kM))]);

for i = 1:4
    if control(i) <= 0
        control(i) = 0;
    elseif control(i) >= 900^2
        control(i) = 900^2;
        control(i) = sqrt(control(i));
    else
        control(i) = sqrt(control(i));
    end
end

end