% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear all;
close;

%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 u1 real

m1 = 20;              % giunto 1 mass
m2 = 20;              % giunto 2 mass
l1 = 0.5;             % lunghezza primo giunto
lc1 = 0.25;           % centro di massa giunto 1
l2 = 0.5;
lc2 = 0.25;
I1 = 0.05;
I2 = 0.05;
G = 9.81;            % constante di accelerazione gravitazionale
f1 = 0.5;
f2 = 0.5;

time = 500;
dt = 0.01;
q1_0 = pi-pi/18;
q2_0 = pi/18;
dq1_0 = 0;
dq2_0 = 0;


%% model without uncertainties

a1 = I1+m1*lc1^2+I2+m2*(l1^2+lc2^2);
a2 = m2*l1*lc2;
a3 = I2+m2*lc2^2;
a4 = G*(m1*lc1+m2*l1);
a5 = G*m2*lc2;

COEFF = [a1 a2 a3 a4 a5 f1 f2];

M = [a1 + 2*a2*cos(q2), a3+a2*cos(q2);
         a3+a2*cos(q2),          a3];
 

c = [a2*sin(q2)*dq2*(dq2+2*dq1);
         a2*sin(q2)*dq1^2];

e = [a4*sin(q1)+a5*sin(q1+q2);
            a5*sin(q1+q2)];
        
n = c+e;

F = diag([f1,f2]);   %friction matrix

f =[dq1;
    dq2;
    -inv(M)*(n+F*[dq1;dq2])
    ];

g =[0,0;
    0,0;
    inv(M)];

dx=f+g*[u1;0];
%% LQR INIT

Q1E = pi;

q1e = Q1E;
q2e = pi - Q1E;

Q_DES = [q1e; q2e];

n = [1 0]';

Q = diag([10 10 1 1]);
R = 1;

tau_eq = a4*sin(q1e)+a5*sin(q1e+q2e);

Me = [a1 + 2*a2*cos(q2e), a3+a2*cos(q2e);
         a3+a2*cos(q2e),          a3];

He = [a4*cos(q1e)+a5*cos(q1e+q2e),  a5*cos(q1e+q2e);
           a5*cos(q1e+q2e),         a5*cos(q1e+q2e)];
       
A = [zeros(2), eye(2);
     -inv(Me)*He,   -inv(Me)*F];
 
b = [0; 0; Me\n];

[K,s,e] = lqrd(A,b,Q,R,0.005);

K_p = -K(1:2);
K_d = -K(3:4);

%% exec
x = zeros(4,time);
dstate = zeros(4,time);
u = zeros(2,time);
dstate(:,1) =subs(dx,[q1,q2,dq1,dq2,u1],[q1_0,q2_0,dq1_0,dq2_0,0]);
x(:,1) = [q1_0,q2_0,dq1_0,dq2_0]';

%%%%%%%%%%%%%%%%%%%%
fig = figure();
j1x = l1*sin(x(1,1));
j1y = -l1*cos(x(1,1));
j2x = j1x + l2*sin(x(2,1)+x(1,1));
j2y = j1y - l2*cos(x(2,1)+x(1,1));

hold on
joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
tip = plot(j2x,j2y,'*','MarkerFaceColor','red','Color','red');
xlim([-1 1]);
ylim([-1,1]); 
hold off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use_cbf = 1;     %flag to use CBF
use_disturb = 1; %flag to introduce disturbance

for i = 1:time
    tau = K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i)+ tau_eq;
    if use_disturb == 1 && i > 200 && i <=220
        disturb = -20;
    else 
        disturb = 0;
    end
    
    if use_cbf == 1
        u(:,i) =  CBFcontroller(COEFF,x(:,i),Q1E,tau);
    else
        u(1,i) = tau;
    end

    dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1],[x(:,i)',u(1,i)+disturb]);
    x(:,i+1) = x(:,i) + dt*dstate(:,i);
    %disp("x(:,i+1) = "+x(:,i)+" + "+dt+"*"+dstate(:,i)+" :");  
    %disp(x(:,i+1));
    %t(i+1) = i*dt;
 
    figure(fig);
    delete(joint1);
    delete(joint2);
    delete(tip);
    j1x = l1*sin(x(1,i));
    j1y = -l1*cos(x(1,i));
    j2x = j1x + l2*sin(x(2,i)+x(1,i));
    j2y = j1y - l2*cos(x(2,i)+x(1,i));

    hold on
    joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
    joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
    tip = plot(j2x,j2y,'*','MarkerFaceColor','red','Color','red');
    xlim([-1 1]);
    ylim([-1,1]);
    hold off
    drawnow();
end

close(fig);

%% plot
%PenduPlot(x,time,l1,l2);

% show the state of the robot

N = (1 : time) * dt;

figure(1)
plot(N,x(1,1:time),N,x(2,1:time));
xlabel('time');
ylabel('q')
legend('q1','q2');
title('q angle')

figure(2)
plot(N,x(3,1:time),N,x(4,1:time));
xlabel('time');
ylabel('dq');
legend('dq1','dq2');
title('dq velocity');

figure(3)
plot(x(1,1:time),x(3,1:time));
xlabel('q1');
ylabel('dq1')
title('q1 and dq1')

figure(4)
plot(x(2,1:time),x(4,1:time));
xlabel('q2');
ylabel('dq2')
title('q2 and dq2')

figure(5)
plot(N,u(1,1:time));
xlabel('time');
ylabel('tau1');
title('input torque');

function u = CBFcontroller(coeff,x,q1e,u_des)

a1 = coeff(1);
a2 = coeff(2);
a3 = coeff(3);
a4 = coeff(4);
a5 = coeff(5);
f1 = coeff(6);
f2 = coeff(7);

q1  = x(1);
q2  = x(2);
dq1 = x(3);
dq2 = x(4);

M = [a1 + 2*a2*cos(q2), a3+a2*cos(q2);
         a3+a2*cos(q2),          a3];

c = [a2*sin(q2)*dq2*(dq2+2*dq1);
         a2*sin(q2)*dq1^2];

e = [a4*sin(q1)+a5*sin(q1+q2);
            a5*sin(q1+q2)];
        
n = c+e;

F = diag([f1,f2]);   %friction matrix

f = [dq1;
    dq2;
    -inv(M)*(n+F*[dq1;dq2])
    ];

g = [0,0;
    0,0;
    inv(M)];

%cbf
q2e = pi - q1e;

Q2E_MAX = pi/18;
c = 1;
alpha = 0.1;
cbf = 0.5*(Q2E_MAX^2-(q2-q2e)^2-c*dq2^2);
grad = [0, q2e - q2, 0, -c*dq2];
% cbf = 0.5*(Q2E_MAX^2-(q2+q1-q1e-q2e)^2-c*(dq1^2+dq2^2));
% grad = [q1e - q2 - q1 + q2e, q1e - q2 - q1 + q2e, -c*dq1, -c*dq2];

H = eye(2);
f_qp = -[u_des; 0];
A = -grad*g;
b = grad*f+alpha*cbf;
options = optimoptions('quadprog','Algorithm','active-set','Display','off');

u = quadprog(H,f_qp,A,b,[],[],[],[],[u_des,0]',options);
end