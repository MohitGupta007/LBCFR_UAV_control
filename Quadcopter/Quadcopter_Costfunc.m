%% ============================================================ %%
%%  Quadcopter Cost function
%% ============================================================ %%

function q = Quadcopter_Costfunc(state,control,du,nu,R,xd)

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

Critical_cost = 0;

if z < 0
    Critical_cost = Critical_cost + 1000;
elseif (phi >= 40*pi/180) || (phi <= -40*pi/180)
    Critical_cost = Critical_cost + 1000;
elseif (theta >= 40*pi/180) || (theta <= -40*pi/180)
    Critical_cost = Critical_cost + 1000;
%{    
elseif abs(vx) >= 3
    Critical_cost = Critical_cost + 100;
elseif abs(vy) >= 3
    Critical_cost = Critical_cost + 100;
elseif abs(vz) >= 3
    Critical_cost = Critical_cost + 100;
elseif abs(p) >= 5*pi/180
    Critical_cost = Critical_cost + 100;
elseif abs(q) >= 5*pi/180
    Critical_cost = Critical_cost + 100;
elseif abs(r) >= 5*pi/180
    Critical_cost = Critical_cost + 100;
%}
end

q = 100*(x-xd(1,1))^2 + 100*(y-xd(1,2))^2 + 100*(z-xd(1,3))^2 + 10*(vx^2+vy^2+vz^2) + 0.1*180^2/(pi^2)*(theta^2+phi^2+psi^2)...
    +(1-nu^(-1))/2*du'*R*du + (control-620.6108)'*R*du + 1/2*(control-610.6108)'*R*(control-610.6108) ...
    +Critical_cost;

end