clear; close; clc;
% ROS Setup
M1 = 1; %Kg
M2 = 1; %Kg
L1 = 1; %m
L2 = 1; %m
r1 = 0.45; %m
r2 = 0.45; %m
I1 = 0.084; %Kg.m2
I2 = 0.084; %Kg.m2
g = 9.81; %m/s2

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

lambda = [-3 -3 -4 -4];

K = place(A,B,lambda);
Kp = K(:,1:2);              %[12 0; 0 12]
Kd = K(:,3:4);
O = [0 0; 0 0];
Acl = [O eye(2); -Kp, -Kd] ;
Q = eye(4).*(0.90);
P = lyap(Acl',Q);
Gamma = eye(5).*50;
alpha0 = [M2*L1^2 + M1*r1^2 + M2*r2^2 + I1 + I2;
            M2*L1*r2;
            M2*r2^2 + I2;
            M1*r1 + M2*L1;
            M2*r2].*0.75; 
alpha = alpha0;

rosinit;
j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);
tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
%req.JointVelocities = [0,0];
resp = call(client,req,'Timeout',3);

i = 1;

tic;
t = 0;
t_prev = 0;
theta1_dot_prev = 0;
theta2_dot_prev = 0;
while(t < 10)
t = toc;
% read the joint states
jointData = receive(JointStates);
% inspect the "jointData" variable in MATLAB to get familiar with its structure
% design your state feedback controller in the following
if t >= 8.0
    theta1 = wrapToPi(jointData.Position(1));
    theta2 = wrapToPi(jointData.Position(2));
else
    theta1 = wrapTo2Pi(jointData.Position(1));
    theta2 = wrapTo2Pi(jointData.Position(2));
end
theta1_dot = jointData.Velocity(1);
theta2_dot = jointData.Velocity(2);
a1 = alpha(1);
a2 = alpha(2);
a3 = alpha(3);
a4 = alpha(4);
a5 = alpha(5);

X = [theta1;theta2;theta1_dot;theta2_dot];

theta1_desired = (pi*t^3)/500 - (3*pi*t^2)/100 - t/18014398509481984 + pi;
theta2_desired = (pi*t^3)/1000 - (3*pi*t^2)/200 - t/36028797018963968 + pi/2;

theta1_dot_desired = (3*pi*t^2)/500 - (3*pi*t)/50 - 1/18014398509481984;
theta2_dot_desired = (3*pi*t^2)/1000 - (3*pi*t)/100 - 1/36028797018963968;

theta1_ddot_desired = (3*pi*t)/250 - (3*pi)/50;
theta2_ddot_desired = (3*pi*t)/500 - (3*pi)/100;

feed_foward_input = [theta1_ddot_desired; theta2_ddot_desired];
e = [theta1 - theta1_desired; theta2 - theta2_desired];
e_dot = [theta1_dot - theta1_dot_desired; theta2_dot - theta2_dot_desired];
G = [e; e_dot];

V = feed_foward_input - Kp*e - Kd*e_dot;
Mmat_hat = [a1 + 2*a2*cos(theta2), a3 + a2*cos(theta2);a3 + a2*cos(theta2),                  a3];
Cmat_hat = [-a2*theta2_dot*sin(theta2)*(2*theta1_dot + theta2_dot);a2*theta1_dot^2*sin(theta2)];
Gmat_hat = [- a4*g*sin(theta1) - a5*g*sin(theta1 + theta2); -a5*g*sin(theta1 + theta2)];

U = Mmat_hat*V + Cmat_hat + Gmat_hat;
tau1.Data = U(1);
tau2.Data = U(2);
send(j1_effort,tau1);
send(j2_effort,tau2);

% you can sample data here to be plotted at the end
X1(i) = theta1;
X2(i) = theta2;
X3(i) = theta1_dot;
X4(i) = theta2_dot;
TAU1(i) = U(1);
TAU2(i) = U(2);

Tau1 = U(1);
Tau2 = U(2);

%theta2_ddot = -(I2*Tau1 - I1*Tau2 - I2*Tau2 - L1^2*M2*Tau2 - M1*r1^2*Tau2 + M2*r2^2*Tau1 - M2*r2^2*Tau2 - L1^2*M2^2*g*r2*sin(theta1 + theta2) + L1*M2^2*g*r2^2*sin(theta1) - I1*M2*g*r2*sin(theta1 + theta2) + I2*L1*M2*g*sin(theta1) + I2*M1*g*r1*sin(theta1) + L1*M2*r2*Tau1*cos(theta2) - 2*L1*M2*r2*Tau2*cos(theta2) + L1*M2^2*r2^3*theta1_dot^2*sin(theta2) + L1^3*M2^2*r2*theta1_dot^2*sin(theta2) + L1*M2^2*r2^3*theta2_dot^2*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot^2*cos(theta2)*sin(theta2) + L1^2*M2^2*r2^2*theta2_dot^2*cos(theta2)*sin(theta2) - L1*M2^2*g*r2^2*sin(theta1 + theta2)*cos(theta2) + L1^2*M2^2*g*r2*cos(theta2)*sin(theta1) - M1*M2*g*r1^2*r2*sin(theta1 + theta2) + I1*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta1_dot^2*sin(theta2) + I2*L1*M2*r2*theta2_dot^2*sin(theta2) + M1*M2*g*r1*r2^2*sin(theta1) + 2*L1*M2^2*r2^3*theta1_dot*theta2_dot*sin(theta2) + 2*L1^2*M2^2*r2^2*theta1_dot*theta2_dot*cos(theta2)*sin(theta2) + L1*M1*M2*r1^2*r2*theta1_dot^2*sin(theta2) + 2*I2*L1*M2*r2*theta1_dot*theta2_dot*sin(theta2) + L1*M1*M2*g*r1*r2*cos(theta2)*sin(theta1))/(- L1^2*M2^2*r2^2*cos(theta2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
TIMEEE = t-t_prev;
%theta1_ddot = (theta1_dot - theta1_dot_prev)/(t-t_prev);
%theta2_ddot = (theta2_dot - theta2_dot_prev)/(t-t_prev);
%theta1_dot_prev = theta1_dot;
%theta2_dot_prev = theta2_dot;

Theta_DDot = Mmat_hat\(U - Gmat_hat - Cmat_hat);
theta1_ddot = Theta_DDot(1);
theta2_ddot = Theta_DDot(2);
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
alpha = alpha + alpha_tilda*(t-t_prev);
A1(i) = alpha(1);
A2(i) = alpha(2);
A3(i) = alpha(3);
A4(i) = alpha(4);
A5(i) = alpha(5);

time(i) = t;
t_prev = t;
i = i+1;

end

tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);
% disconnect from roscore
rosshutdown;

Theta1_desired = (pi*time.^3)/500 - (3*pi*time.^2)/100 - time/18014398509481984 + pi;
Theta2_desired = (pi*time.^3)/1000 - (3*pi*time.^2)/200 - time/36028797018963968 + pi/2;

Theta1_dot_desired = (3*pi*time.^2)/500 - (3*pi*time)/50 - 1/18014398509481984;
Theta2_dot_desired = (3*pi*time.^2)/1000 - (3*pi*time)/100 - 1/36028797018963968;


%% plots
figure
hold on
subplot(2,2,1)
plot(time,X1)
hold on
plot(time, Theta1_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad')
title('theta1')

subplot(2,2,2)
plot(time,X2)
hold on
plot(time, Theta2_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad')
title('theta2')

subplot(2,2,3)
plot(time,X3)
hold on
plot(time, Theta1_dot_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad/s')
title('theta1-dot')

subplot(2,2,4)
plot(time,X4)
hold on
plot(time, Theta2_dot_desired,'Color','red','LineStyle','--')
xlabel('Time step')
ylabel('rad/s')
title('theta2-dot')
hold off

figure
hold on
subplot(2,1,1)
plot(time,TAU1)
xlabel('Time step')
ylabel('Nm')
title('tau1')

subplot(2,1,2)
plot(time,TAU2)
xlabel('Time step')
ylabel('Nm')
title('tau2')
hold off

figure
hold on
subplot(2,3,1)
plot(time,A1)
xlabel('Time step')
ylabel('a1')
title('a1')

subplot(2,3,2)
plot(time,A2)
xlabel('Time step')
ylabel('a2')
title('a2')

subplot(2,3,3)
plot(time,A3)
xlabel('Time step')
ylabel('a3')
title('a3')

subplot(2,3,4)
plot(time,A4)
xlabel('Time step')
ylabel('a4')
title('a4')

subplot(2,3,5)
plot(time,A5)
xlabel('Time step')
ylabel('a5')
title('a5')
