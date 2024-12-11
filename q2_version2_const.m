x0 = [pi/3;pi/3;0;0;0;0;0;0;0];
Ts = [0 30];
lamda= 2*eye(2);
[T,X] =ode45(@(t,x) model2dof(t,x,lamda),Ts,x0);
getplots(T,X,lamda);


function  xdot = model2dof(t,x,lamda)
%Specifications
m1 = 6;
m2 = 4;
l1 = 0.5;
l2 = 0.4;
g = 9.81;
lc1 = 0.2;
lc2 = 0.1;
I1 = 0.43;
I2 = 0.05;
ml = 0.5;
qd = [pi/2;-pi/3];
qddot = [0;0];
qddotdot = [0;0];
q =[x(1);x(2)];
qdot  = [x(3);x(4)];
%definition of error
qbar = q - qd;
qbardot = qdot - qddot;
%new reference trajectory
qr =  qddot - lamda*qbar;
qrdot = qddotdot - lamda*qbardot;
%definition of s
s = qbardot + lamda*qbar;
%H matrix
h11 = m1*lc1^2 + m2*((lc2)^2+(l1)^2+2*l1*lc2*cos(x(2))) + ml*(l2^2+l1^2+2*l1*l2*cos(x(2)))+I1+I2;
h12 = m2*lc2*(lc2+l1*cos(x(2)))+ ml*l2*(l2+l1*cos(x(2)))+I2;
h22 = (lc2^2)*m2+(l2^2)*ml+I2;
H = [h11 h12;h12 h22];
%C matrix 
c11 = -l1*(m2*lc2+ml*l2)*sin(x(2))*x(4);
c12 = -l1*(m2*lc2+ml*l2)*sin(x(2))*(x(3)+x(4));
c21 = l1*(m2*lc2+ml*l2)*sin(x(2))*x(3);
C = [c11 c12;c21 0];
% g matrix
g1 = (m2*lc2 + ml*l2)*g*cos(x(1)+x(2))+(m2*l1+ml*l1+m1*lc1)*g*cos(x(1));
g2 = (m2*lc2 + ml*l2)*g*cos(x(1)+x(2));
G = [g1;g2];
%getting the vector of parameters
p = x(5:9);
if(t>29.99)
    p
end
%H estimation
h11e = p(1) + 2*p(2)*cos(x(2));
h12e = p(3) + p(2)*cos(x(2));
h22e = p(3);
He = [h11e h12e;h12e h22e];
%C estimation 
c11e = -p(2)*sin(x(2))*x(4);
c12e = -p(2)*sin(x(2))*(x(4)+x(3));
c21e = p(2)*sin(x(2))*x(3);
Ce =[c11e c12e;c21e 0];
%G estimation
g1e = p(4)*g*cos(x(1)+x(2)) + p(5)*g*cos(x(1));
g2e = p(4)*g*cos(x(1)+x(2));
Ge = [g1e;g2e];
%Y matrix (the transpose)
Yt = [ qrdot(1) 0;
    2*qrdot(1)*cos(x(2))+qrdot(2)*cos(x(2))-qr(1)*sin(x(2))*(qdot(1)+qdot(2)) qrdot(1)*cos(x(2))+sin(x(2))*qr(1)*qdot(1);
    qrdot(2) qrdot(1)+qrdot(2);
    g*cos(x(1)+x(2)) g*cos(x(1)+x(2));
    g*cos(x(1)) 0
   ];
%definition of K
K = 2*eye(2);
%definition the u
u = He*qrdot + Ce*qr + Ge - K*s ;
%cast(norm(s),'single')
qdd = inv(H)*(u-G)- inv(H)*C*qd;
Ginv = 5*eye(5);
dparam = -Ginv*Yt*s;
dparamt = dparam';
t
%Definition of state space of the system
xdottemp = [x(3);x(4);qdd(1);qdd(2)];
xdot=[xdottemp;dparam];
end
function getplots(T,X,lamda) 
    qd = [pi/2;-pi/3];
    qd1 = @(t)(pi/2);
    qd2 = @(t)(-pi/3);
    qddot = [0;0];
    qddot1 = @(t)(0);
    qddot2 = @(t)(0);
    qddotdot = [0;0];
    qddotdot1 = @(t)(0);
    qddotdot2 = @(t)(0);
    q1= X(:,1);
    q2 = X(:,2);
    qdot1 = X(:,3);
    qdot2 = X(:,4);
   
%definition of error
    qbar1 = q1 - qd1(T);
    qbardot1 =qdot1 - qddot1(T);
    qbar2 =q2 - qd2(T);
    qbardot2 =qdot2- qddot2(T);
%new reference trajectory
    lamda1 = lamda(1,1);
    lamda2 = lamda(2,2);
    qr1 =  qddot1(T) - lamda1*qbar1;
    qrdot1 = qddotdot1(T) - lamda1*qbardot1;
    qr2 =  qddot2(T) - lamda2*qbar2;
    qrdot2 = qddotdot2(T) - lamda2*qbardot2;
    e1 = X(:,1) -qd1(T);
    e2 = X(:,2) -qd2(T) ; 
    s1 =  qbardot1 + lamda1*qbar1;
    s2 =  qbardot2 + lamda2*qbar2; 
    
    figure(1);
    plot(T,X(:,1));
    title('Angular Position 1')
    xlabel('Time in sec')
    ylabel('Angle in rad')
    grid;
    hold;
    figure(2);
    plot(T,X(:,2));
    title('Angular Position 2')
    xlabel('Time in sec')
    ylabel('Angle in rad')
    grid;
    figure(3);
    plot(T,X(:,3));
    title('Angular Velocity 1')
    xlabel('Time in sec')
    ylabel('Anglular velocity in rad/s')
    grid;
    hold ;
    figure(4);
    plot(T,X(:,4));
    title('Angular Velocity 2')
    xlabel('Time in sec')
    ylabel('Anglular velocity in rad/s')
    grid;
    hold;
    figure(5);
    plot(T,qbar1);
    title('Error for anglular position 1,qbar1')
    xlabel('Time in sec')
    ylabel(' Error of anglular position in rad')
    grid;
    hold;
    figure(6);
    plot(T,qbar2);
    title('Error for anglular position 2,qbar2')
    xlabel('Time in sec')
    ylabel(' Error of anglular position in rad')
    grid;
    hold;
    figure(7);
    plot(T,s1);
    title('S1 of addaptive control')
    xlabel('Time in sec')
    grid;
    hold;
    figure(8);
    plot(T,s2);
    title('S2 of addaptive control')
    xlabel('Time in sec')
    grid;
    hold;
    figure(9);
    plot(T,qr1);
    title('Qr1')
    xlabel('Time in sec')
    ylabel(' Angle position in rad')
    grid;
    hold;
    figure(10);
    plot(T,qr2);
    title('Qr2')
    xlabel('Time in sec')
    ylabel('Angle position in rad')
    grid;
    hold;
    

end


