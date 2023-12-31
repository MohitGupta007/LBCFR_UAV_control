function next_state = get_new_states(state,control,dt,model)

    g  = 9.81;
    m  = 0.468;
    L  = 0.225;
    km = 20;
    kF = 2.980*10^(-6);
    kM = 1.140*10^(-7);
    IM = 3.357*10^(-5);
    Ixx = 4.856*10^(-3);
    Iyy = 4.856*10^(-3);
    Izz = 8.801*10^(-3);
    Ax = 0.25;
    Ay = 0.25;
    Az = 0.25;

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
    sigma1  = state(13);
    sigma2  = state(14);
    sigma3  = state(15);
    sigma4  = state(16);
    
    Position     = [x;y;z];
    Velocity     = [vx;vy;vz];
    Angle        = [phi;theta;psi];
    AngularVel   = [p;q;r];
    Sigma        = [sigma1;sigma2;sigma3;sigma4];

    F = kF*[sigma1^2; sigma2^2; sigma3^2; sigma4^2];
    M = kM*[sigma1^2; sigma2^2; sigma3^2; sigma4^2];
    Tau = [L*(-F(2)+F(4)); L*(-F(1)+F(3)); -M(1)+M(2)-M(3)+M(4)];
    T = sum(F);
    Gamma = sigma1-sigma2+sigma3-sigma4;
    
    x_dot = model.predict([x,y,z,vx,vy,vz,phi,theta,psi,p,q,r,T,Tau',Gamma])';
    
    next_state        = zeros(12,1);
    next_state(1:3)   = Position   + x_dot(1:3)*dt;
    next_state(4:6)   = Velocity   + x_dot(4:6)*dt;
    next_state(7:9)   = Angle      + x_dot(7:9)*dt;
    next_state(10:12) = AngularVel + x_dot(10:12)*dt;
    next_state(13:16) = Sigma      + km*(control-Sigma)*dt;
end