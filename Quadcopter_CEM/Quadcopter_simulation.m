%% ============================================================ %%
%%  Quadcopter Simulation
%% ============================================================ %%

clear all;

model = importTensorFlowNetwork('model_0405_8_lr');

%%

% Quadcopter state
state       = zeros(16,1);
state_diff  = zeros(16,1);
state(3)    = 10;
state(7)    = 0*pi/180;
state(8)    = 0*pi/180;
state(13:16)= 620.6108;

N = 3000;
state_log = zeros(16,N);
state_diff_log = zeros(16,N);
control_log = zeros(4,N);
%T_and_Tau_log = zers(4,N);
control = 620.6108*ones(4,1);

params.plant_dt = 0.01;
params.control_dt = 0.05;
params.Hp = 15;
params.K = 200;
params.num_iterations = 4;

params.sigma = 100.*ones(4,params.Hp);

% circular path
t = 0:params.plant_dt:(N-1)*params.plant_dt;
r = 2;
omega = pi/10;
xc = r*cos(omega*t);
yc = r*sin(omega*t);
zc = 11*ones(size(t));
xd = [xc',yc',zc'];

% % single waypoint path
% t = 0:params.plant_dt:(N-1)*params.plant_dt;
% xc = 3*ones(size(t));
% yc = 3*ones(size(t));
% zc = 13*ones(size(t));
% xd = [xc',yc',zc'];

traj_err = [];

tic
for i = 1:N
    state_log(:,i)  = state;
    state_diff_log(:,i)  = state_diff;
    
    if((i+params.Hp-1)<=N)
        params.xd = xd(i:i+params.Hp-1,:);
    else
        params.xd(1:params.Hp,:) = repmat(xd(end,:),params.Hp,1);
        params.xd(1:min(i+params.Hp-1,size(xd,1))-i+1,:) = xd(i:min(i+params.Hp-1,size(xd,1)),:);
    end
%     tic
    if rem(i-1,params.control_dt/params.plant_dt) == 0
        i
        if i == 1
            [control,sigma] = Quadcopter_CE(state,model,params);
        else
            [control,sigma] = Quadcopter_CE(state,model,params,control);
        end
    end
%     toc
%     if rem(i-1,2) == 0
%         [control, T_and_Tau] = Quadcopter_PID_controller(state,state_diff,control);
%     end

    state = Quadcopter_Dynamics(state,control(:,1),params.plant_dt);
    err = norm(params.xd(1,:)-state(1:3)')
    traj_err = [traj_err norm(params.xd(1,:)-state(1:3)')];
    state_diff = (state-state_log(:,i))/params.plant_dt;
    control_log(:,i) = control(:,1);
    %T_and_Tau_log(:,i) = T_and_Tau;
    
end
toc
%%
close all;
plant_dt = 1;

Quadcopter_Animator(state_log);
hold on
plot3(state_log(1,:),state_log(2,:),state_log(3,:),'b')
% scatter3(xd(:,1),xd(:,2),xd(:,3),20,'r','filled','MarkerEdgeColor','k')
plot3(xd(:,1),xd(:,2),xd(:,3),'r')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
grid on
title('Quadcopter Simulation')
hold off

figure
tiledlayout(3,1)
nexttile
plot(t.*plant_dt,xd(:,1))
hold on
plot(t.*plant_dt,state_log(1,:))
grid on
xlabel('time (s)')
ylabel('X (m)')
title('Position Tracking')
hold off
nexttile
plot(t.*plant_dt,xd(:,2))
hold on
plot(t.*plant_dt,state_log(2,:))
grid on
xlabel('time (s)')
ylabel('Y (m)')
hold off
nexttile
plot(t.*plant_dt,xd(:,3))
hold on
plot(t.*plant_dt,state_log(3,:))
grid on
xlabel('time (s)')
ylabel('Z (m)')
hold off

figure
tiledlayout(6,1)
nexttile
plot(t.*plant_dt,state_log(4,:))
nexttile
plot(t.*plant_dt,state_log(5,:))
nexttile
plot(t.*plant_dt,state_log(6,:))
nexttile
plot(t.*plant_dt,state_log(10,:))
nexttile
plot(t.*plant_dt,state_log(11,:))
nexttile
plot(t.*plant_dt,state_log(12,:))

figure
tiledlayout(4,1)
nexttile
plot(t.*plant_dt,control_log(1,:)')
grid on
title('Control signals')
xlabel('time (s)')
ylabel('sigma1 (rpm)')
nexttile
plot(t.*plant_dt,control_log(2,:))
grid on
xlabel('time (s)')
ylabel('sigma2 (rpm)')
nexttile
plot(t.*plant_dt,control_log(3,:))
grid on
xlabel('time (s)')
ylabel('sigma3 (rpm)')
nexttile
plot(t.*plant_dt,control_log(4,:))
grid on
xlabel('time (s)')
ylabel('sigma4 (rpm)')

figure
plot(t.*plant_dt,vecnorm(control_log,2,1)')
grid on
title('Control effort')
xlabel('time (s)')
ylabel('effort')

figure
plot(t.*plant_dt,traj_err)
grid on
xlabel('time (s)')
ylabel('Error (m)')
title('Trajectory error')
%% ============================================================ %%
%% NewFile and plotdata(png) creating (Taking account of the date and time)            
%% ============================================================ %%
%{
% Newfile Name
FileName = strrep(strrep(strcat('QuadcoptermatFile_',datestr(datetime('now'))),':','_'),' ','_');

% present date and time
date = strrep(strrep(strcat('_',datestr(datetime('now')),'.png'),':','_'),' ','_');
mat = strrep(strrep(strcat('_',datestr(datetime('now')),'.mat'),':','_'),' ','_');

% creating Newfile
mkdir(FileName);

% creating pngfile from plotfigure
save([FileName,'/result',mat]);
%}