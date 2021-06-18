% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear all;
close;

%% configuration

use_cbf = true;
use_lqr = false;
disturbance = false;
show_animation = true;
show_plots = true;

%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 u1 real
global Q_DES  Q2_MAX
global a1 a2 a3 a4 a5 f1 f2 

m1 = 2;              % giunto 1 mass
m2 = 2;              % giunto 2 mass
l1 = 0.5;             % lunghezza primo giunto
lc1 = 0.25;           % centro di massa giunto 1
l2 = 0.5;
lc2 = 0.25;
I1 = 0.05;
I2 = 0.05;
G = 9.81;            % constante di accelerazione gravitazionale
f1 = 0.5;
f2 = 0.5;

time = 1000;
dt = 0.01;
q1_0 = pi-0.05;
q2_0 = 0;
dq1_0 = 0;
dq2_0 = 0;

Q2_MAX = pi/12;

%% model without uncertainties
a1 = I1+m1*lc1^2+I2+m2*(l1^2+lc2^2);
a2 = m2*l1*lc2;
a3 = I2+m2*lc2^2;
a4 = G*(m1*lc1+m2*l1);
a5 = G*m2*lc2;


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

    Q2E = pi - Q1E;

    Q_DES = [Q1E; Q2E];

    n = [1 0]';

    Q = diag([10 100 1 1]);
    R = 1;

    tau_eq = a4*sin(Q1E)+a5*sin(Q1E+Q2E);

    Me = [a1 + 2*a2*cos(Q2E), a3+a2*cos(Q2E);
             a3+a2*cos(Q2E),          a3];

    He = [a4*cos(Q1E)+a5*cos(Q1E+Q2E),  a5*cos(Q1E+Q2E);
               a5*cos(Q1E+Q2E),         a5*cos(Q1E+Q2E)];

    A = [zeros(2), eye(2);
         -inv(Me)*He,   -inv(Me)*F];

    b = [0; 0; Me\n];
if(use_lqr)
    
    [K,s,e] = lqr(A,b,Q,R);

    K_p = -K(1:2);
    K_d = -K(3:4);
else
    K_p = zeros(1,2);
    K_d = zeros(1,2);
end

%% exec


x = zeros(4,time);
dstate = zeros(4,time);

u = zeros(2,time);              %the final input
tau = zeros(1,time);            %LQR output

h = zeros(1,time);
q_des = zeros(2,time);
q_des(1,:)=pi; 
dstate(:,1) = subs(dx,[q1,q2,dq1,dq2,u1],[q1_0,q2_0,dq1_0,dq2_0,0]);
x(:,1) = [q1_0,q2_0,dq1_0,dq2_0]';

if(show_animation == 1)
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
    drawnow;
end

for i = 1:time
    
    Q_DES(2) = Q2E + (pi/144)*sin(i*pi/250);

    tau(i) = K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i) + tau_eq;

    if use_cbf == 1
        [h(i), u(1,i)] =  CBFcontroller(x(:,i),tau(i));
%         mpc_x = [x(:,i)];
%         mpc_tau = [];
%         for idx = 1:20
%             mpc_tau = [mpc_tau; K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i) + tau_eq]; 
%             mpc_dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1],[x(:,i)',u(1,i)]);
%             mpc_x  =[mpc, mpc_x(:,i+idx) + dt*mpc_dstate(:,i)];
%         end
%     
    else
        u(1,i) = tau(i);
    end

    dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1],[x(:,i)',u(1,i)]);
    x(:,i+1) = x(:,i) + dt*dstate(:,i);
    %disp("x(:,i+1) = "+x(:,i)+" + "+dt+"*"+dstate(:,i)+" :");  
    %disp(x(:,i+1));
    %t(i+1) = i*dt;
 
    if show_animation == 1
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
end

%% plot
% show the state of the robot
if show_plots
    N = (1 : time) * dt;

    figure();
    subplot(2,2,1);
    plot(N,x(2,1:time),N,q_des(2,1:time))
    yline(Q2_MAX);
    yline(-Q2_MAX);
    xlabel('s')
    ylabel('q2')
    legend('q2','q2_des');
    title('q2 angle')

    subplot(2,2,2);
    plot(N,x(1,1:time),N,q_des(1,1:time))
    yline(Q2_MAX);
    yline(-Q2_MAX);
    xlabel('s')
    ylabel('q1')
    legend('q1','q1_des');
    title('q1 angle')

    subplot(2,2,3)
    plot(N,x(3,1:time),N,x(4,1:time));
    xlabel('time');
    ylabel('dq');
    legend('dq1','dq2');
    title('dq velocity');

    subplot(2,2,4)
    plot(N,u(1,1:time));
    xlabel('t');
    ylabel('dq1')
    title('q1 and dq1')


    if(use_cbf == 1)
        figure;
        subplot(1,2,1)
        plot(N,u(1,:),N,tau)
        xlabel('time')
        ylabel('input torque')
        legend('u','LQR')
        title('Final control input vs LQR component')

        subplot(1,2,2)
        plot(N,h)
        xlabel('time')
        ylabel('h')
        title('Control Barrier function')
    else
        figure(5)
        plot(N,u(1,:))
        xlabel('time')
        ylabel('tau1')
        title('LQR input')
    end
end
%% CBF
function [h, u] = CBFcontroller(x,u_des)
global a1 a2 a3 a4 a5 f1 f2 Q_DES Q2_MAX
syms q1 q2 dq1 dq2;

q1  = x(1);
q2  = x(2);
dq1 = x(3);
dq2 = x(4);

Q1E = Q_DES(1);
Q2E = Q_DES(2);

n = [1 0]';

M = [a1 + 2*a2*cos(q2), a3+a2*cos(q2);
         a3+a2*cos(q2),          a3];

c = [a2*sin(q2)*dq2*(dq2+2*dq1);
         a2*sin(q2)*dq1^2];

e = [a4*sin(q1)+a5*sin(q1+q2);
            a5*sin(q1+q2)];

F = diag([f1,f2]);   %friction matrix

f = [dq1;
    dq2;
    -M\(c+e+F*[dq1;dq2])
    ];

g = [0; 0; M\n];

%cbf
%disp(f);

c_param = 1;
alpha = 20;
cbf = pi- q1-q2- pi/12;
grad = [-1, -1, 0, 0];
 
del =pi - pi/10;
cbf = del^2 +q1^2+ 2*q1*q2 +q2^2 - 2*del*(q1+q2);

grad =[2+2*q2-2*del*q2,2+2*q1-2*del*q1,0 ,0 ];


H = 1;
f_qp = -u_des;

A = -grad*g;
b = grad*f+alpha*cbf;

disp("A:"); disp(A);
disp("b:"); disp(b);
options = optimset('display','off');
h = cbf;
u = quadprog(H,f_qp,A,b,[],[],[],[],[],options)
if isempty(u)
    u = 0;
end
end