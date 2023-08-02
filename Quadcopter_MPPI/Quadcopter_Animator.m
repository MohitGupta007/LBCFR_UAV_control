%% ============================================================ %%
%%  Quadcopter Animator
%% ============================================================ %%

function [] = Quadcopter_Animator(state_log)

x     = state_log(1,:);
y     = state_log(2,:);
z     = state_log(3,:);
phi   = state_log(7,:);
theta = state_log(8,:);
psi   = state_log(9,:);

L     = 0.225;
p_d   = 0.15;

for i = 1:length(x)
    if rem(i,5) == 0
        center = [x(i); 
                  y(i);
                  z(i)];
        vertex1_b = [L;
                     0;
                     0];
        vertex2_b = [0;
                    -L;
                    0];
        vertex3_b = [-L;
                      0;
                      0];
        vertex4_b = [0;
                     L;
                     0];

        propeller1_b = [];
        propeller2_b = [];
        propeller3_b = [];
        propeller4_b = [];
    
        for t = 0:360
            subpropeller1 = [L+p_d*cos(t*pi/180);
                             p_d*sin(t*pi/180);
                             0];
            subpropeller2 = [p_d*cos(t*pi/180);
                            -L+p_d*sin(t*pi/180);
                             0]; 
            subpropeller3 = [-L+p_d*cos(t*pi/180);
                             p_d*sin(t*pi/180);
                             0];
            subpropeller4 = [p_d*cos(t*pi/180);
                             L+p_d*sin(t*pi/180);
                             0];         
            propeller1_b = horzcat(propeller1_b,subpropeller1);
            propeller2_b = horzcat(propeller2_b,subpropeller2);
            propeller3_b = horzcat(propeller3_b,subpropeller3);
            propeller4_b = horzcat(propeller4_b,subpropeller4);
        end
    
        R = [cos(psi(i))*cos(theta(i)) cos(psi(i))*sin(theta(i))*sin(phi(i))-sin(psi(i))*cos(phi(i)) cos(psi(i))*sin(theta(i))*cos(phi(i))+sin(psi(i))*sin(phi(i));
             sin(psi(i))*cos(theta(i)) sin(psi(i))*sin(theta(i))*sin(phi(i))+cos(psi(i))*cos(phi(i)) sin(psi(i))*sin(theta(i))*cos(phi(i))-cos(psi(i))*sin(phi(i));
             -sin(theta(i)) cos(theta(i))*sin(phi(i)) cos(theta(i))*cos(phi(i))];

        propeller1_i = R*propeller1_b+center;
        propeller2_i = R*propeller2_b+center;
        propeller3_i = R*propeller3_b+center;
        propeller4_i = R*propeller4_b+center;

        vertex1_i = R*vertex1_b+center;
        vertex2_i = R*vertex2_b+center;
        vertex3_i = R*vertex3_b+center;
        vertex4_i = R*vertex4_b+center;

        arm1x = [center(1),vertex1_i(1)];
        arm1y = [center(2),vertex1_i(2)];
        arm1z = [center(3),vertex1_i(3)];

        arm2x = [center(1),vertex2_i(1)];
        arm2y = [center(2),vertex2_i(2)];
        arm2z = [center(3),vertex2_i(3)];

        arm3x = [center(1),vertex3_i(1)];
        arm3y = [center(2),vertex3_i(2)];
        arm3z = [center(3),vertex3_i(3)];

        arm4x = [center(1),vertex4_i(1)];
        arm4y = [center(2),vertex4_i(2)];
        arm4z = [center(3),vertex4_i(3)];

        scatter3(5,5,15,15,'r','filled','MarkerEdgeColor','k');
        hold on;
        axis equal;
        grid on;
        plot3(arm1x,arm1y,arm1z,'LineWidth',1.5,'Color','b');
        hold on;
        plot3(arm2x,arm2y,arm2z,'LineWidth',1.5,'Color','b');
        hold on;
        plot3(arm3x,arm3y,arm3z,'LineWidth',1.5,'Color','b');
        hold on;
        plot3(arm4x,arm4y,arm4z,'LineWidth',1.5,'Color','b');
        hold on;

        plot3(propeller1_i(1,:),propeller1_i(2,:),propeller1_i(3,:),'LineWidth',1.5,'Color','r');
        hold on;
        plot3(propeller2_i(1,:),propeller2_i(2,:),propeller2_i(3,:),'LineWidth',1.5,'Color','k');
        hold on;
        plot3(propeller3_i(1,:),propeller3_i(2,:),propeller3_i(3,:),'LineWidth',1.5,'Color','k');
        hold on;
        plot3(propeller4_i(1,:),propeller4_i(2,:),propeller4_i(3,:),'LineWidth',1.5,'Color','k');
        hold on;
        xlim([min(x)-L-p_d max(x)+L+p_d]);
        ylim([min(y)-L-p_d max(y)+L+p_d]);
        zlim([min(z)-L-p_d max(z)+L+p_d]);
        hold off;
        drawnow;
    end
end

end

    
