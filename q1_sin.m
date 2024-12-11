x0 = [pi/3;pi/3;0;0];
Ts = [0 20];
lamda1 =10;
lamda2 =10;
lamda = [lamda1 0;0 lamda2]; 
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
%xdesired
xd = [pi/4+(pi/6)*sin(0.2*pi*t);-pi/3+(pi/3)*cos(0.2*pi*t)];
xddot = [pi^2*cos(pi*t/5)/30;-pi^2*sin(pi*t/5)/15];
xddotdot = [-pi^3*sin(pi*t/5)/150;-pi^3*cos(pi*t/5)/75];
%calculation for controller
normHmax = 3.0267;
normCmax = 0.6364*sqrt(x(3)^2+x(4)^2);
normGmax = 1.995*g;
%definition of errors
e = [x(1)-xd(1);x(2)-xd(2)];
edot=[x(3)-xddot(1);x(4)-xddot(2)];
%definition of surfaces for sliding method 
s =lamda*e + edot;
%definition of estimation of matrices 
Hestimation = [1.945+1.1*cos(x(2)) 0.425+0.55*cos(x(2));0.425+0.55*cos(x(2)) 0.425];
Cestimation = 0.55*[-x(4)*sin(x(2)) -(x(3)+x(4))*sin(x(2));x(3)*sin(x(2)) 0];
Gestimation = [1.1*g*cos(x(1)+x(2))+4*g*cos(x(1));1.1*g*cos(x(1)+x(2))];
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
qd = [x(3);x(4)];
r1 = 1*(-normCmax*norm(qd)- normGmax -normHmax*norm(xddotdot) - normHmax*norm(lamda)*norm(edot)) -1;
r2 =  0.5*(-normCmax)-1;
epsilon = 0.05;
t
if norm(s) >epsilon
    u =  Cestimation*qd + Gestimation + Hestimation*xddotdot - Hestimation*lamda*edot - 0.5*Cestimation*s+r1*(s/norm(s))+ r2*s;
else 
     u =  Cestimation*qd + Gestimation + Hestimation*xddotdot - Hestimation*lamda*edot - 0.5*Cestimation*s+r1*(s/epsilon)+ r2*s;
end 
%{ 
%remove coments for maxtorque implementation
maxtorque = 50;
if u(1) >maxtorque
    u(1) =maxtorque;
end
if u(2)>maxtorque
    u(2)=maxtorque;
end
if u(1) <-maxtorque
    u(1) =-maxtorque;
end
if u(2) <-maxtorque
    u(2) = -maxtorque;
end
%}
%u = [u1;u2];
%cast(norm(s),'single')
qdd = inv(H)*(u-G)- inv(H)*C*qd;
%Definition of state space of the system
xdot = [x(3);x(4);qdd(1);qdd(2)];
end
function getplots(T,X,lamda) 
    qd = @(t)[pi/4+(pi/6)*sin(0.2*pi*t);-pi/3+(pi/3)*cos(0.2*pi*t)];
    qd1 = @(t)(pi/4+(pi/6)*sin(0.2*pi*t));
    qd2 = @(t)(-pi/3+(pi/3)*cos(0.2*pi*t));
    qddot1 =@(t)(pi^2*cos(pi*t/5)/30);
    qddot2 = @(t)(-pi^2*sin(pi*t/5)/15);
    e1 = X(:,1) -qd1(T) ;
    e2 = X(:,2) -qd2(T) ;
    e1dot = X(:,3) - qddot1(T);
    e2dot = X(:,4) - qddot2(T);
    lamda1 = lamda(1,1)
    lamda2 = lamda(2,2);
    s1 = e1dot + lamda1*e1;
    s2 = e2dot + lamda2*e2;
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
    hold;
    figure(4);
    plot(T,X(:,4));
    title('Angular Velocity 2')
    xlabel('Time in sec')
    ylabel('Anglular velocity in rad/s')
    grid;
    hold;
    figure(5);
    plot(T,e1);
    title('Error for anglular position 1')
    xlabel('Time in sec')
    ylabel(' Error of anglular position in rad')
    grid;
    hold;
    figure(6);
    plot(T,e2);
    title('Error for anglular position 2')
    xlabel('Time in sec')
    ylabel(' Error of anglular position in rad')
    grid;
    hold;
    figure(7);
    plot(T,s1);
    title('Sliding  surface 1')
    xlabel('Time in sec')
    grid;
    hold;
    figure(8);
    plot(T,s2);
    title('Sliding surface 2')
    xlabel('Time in sec')
    grid;
    hold;
end


