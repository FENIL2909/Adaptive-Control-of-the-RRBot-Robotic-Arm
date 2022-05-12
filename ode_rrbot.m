%% ode_rrbot
function dX = ode_rrbot(t,X)
Adaptive_Control = true;
% M1 = 1; %Kg
% M2 = 1; %Kg
% L1 = 1; %m
% L2 = 1; %m
% r1 = 0.45; %m
% r2 = 0.45; %m
% I1 = 0.084; %Kg.m2
% I2 = 0.084; %Kg.m2
g = 9.81; 

%m/s2

dX = zeros(9,1);
X = num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot, a1,a2,a3,a4,a5] = deal(X{:});

theta1_desired = (pi*t^3)/500 - (3*pi*t^2)/100 - t/18014398509481984 + pi;
theta2_desired = (pi*t^3)/1000 - (3*pi*t^2)/200 - t/36028797018963968 + pi/2;

theta1_dot_desired = (3*pi*t^2)/500 - (3*pi*t)/50 - 1/18014398509481984;
theta2_dot_desired = (3*pi*t^2)/1000 - (3*pi*t)/100 - 1/36028797018963968;

theta1_ddot_desired = (3*pi*t)/250 - (3*pi)/50;
theta2_ddot_desired = (3*pi*t)/500 - (3*pi)/100;

feed_foward_input = [theta1_ddot_desired; theta2_ddot_desired];
e = [theta1 - theta1_desired; theta2 - theta2_desired];
e_dot = [theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];
%% Virtual Control input design

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
lambda = [-3 -3 -4 -4];

K = place(A,B,lambda);
Kp = K(:,1:2);            
Kd = K(:,3:4);
O = [0 0; 0 0];
Acl = [O eye(2); -Kp, -Kd] ;
Q = eye(4).*(0.9);
if Adaptive_Control
    P = lyap(Acl',Q);
else
    P = 0;
end
Gamma = eye(5).*50;
G = [e; e_dot];

V = feed_foward_input - Kp*e - Kd*e_dot;

Mmat_hat = [a1 + 2*a2*cos(theta2), a3 + a2*cos(theta2);a3 + a2*cos(theta2),                  a3];
Cmat_hat = [-a2*theta2_dot*sin(theta2)*(2*theta1_dot + theta2_dot);a2*theta1_dot^2*sin(theta2)];
Gmat_hat = [- a4*g*sin(theta1) - a5*g*sin(theta1 + theta2); -a5*g*sin(theta1 + theta2)];

U = Mmat_hat.*V + Cmat_hat + Gmat_hat;

tau1 = U(1);
tau2 = U(2);

M1 = 1; %Kg
M2 = 1; %Kg
L1 = 1; %m
L2 = 1; %m
r1 = 0.45; %m
r2 = 0.45; %m
I1 = 0.084; %Kg.m2
I2 = 0.084; %Kg.m2
g = 9.81; %m/s2

dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = (I2*tau1 - I2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 + L1*M2^2*g*r2^2*sin(theta1) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) - L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
dX(4) = -(I2*tau1 - I1*tau2 - I2*tau2 - L1^2*M2*tau2 - M1*r1^2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 - L1^2*M2^2*g*r2*sin(theta1 + theta2) + L1*M2^2*g*r2^2*sin(theta1) - I1*M2*g*r2*sin(theta1 + theta2) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) + L1*M2*r2*tau1*cos(theta2) - 2*L1*M2*r2*tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1^3*M2^2*r2*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + L1^2*M2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + L1^2*M2^2*g*r2*cos(theta2)*sin(theta1) - M1*M2*g*r1^2*r2*sin(theta1 + theta2) + I1*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + L1*M1*M2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2) + L1*M1*M2*g*r1*r2*cos(theta2)*sin(theta1))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);

theta1_ddot = dX(3);
theta2_ddot = dX(4);

Y = [theta1_ddot, ...
     cos(theta2)*(2*theta1_ddot + theta2_ddot) - 2*sin(theta2)*theta1_dot*theta2_dot - sin(theta2)*theta2_dot^2, ...
     theta2_ddot, ...
    -sin(theta1)*g, ...
    -sin(theta1 + theta2)*g; ...
     0, ...
     sin(theta2)*theta1_dot^2 + cos(theta2)*theta1_ddot, ...
     theta1_ddot + theta2_ddot, ...
     0, ...
     -sin(theta1+theta2)*g];

Phi = Mmat_hat\Y;
alpha_tilda = -(Gamma\(Phi'*B'*P*G));

dX(5) = alpha_tilda(1);
dX(6) = alpha_tilda(2);
dX(7) = alpha_tilda(3);
dX(8) = alpha_tilda(4);
dX(9) = alpha_tilda(5);
end
