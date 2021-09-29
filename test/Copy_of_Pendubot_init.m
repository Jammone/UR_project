% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

%clc;
clear all;

%% configuration

use_cbf = true;
use_lqr = true;
show_animation = false;
show_plots = false;

%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 u1 m2_unc real
global Q_DES Q1E Q2E Q1_MAX Q2_MAX
global a1 a2 a3 a4 a5 f1 f2 diffmax



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
dt = 0.005;
q1_0 = pi+pi/18;
q2_0 = pi/18;
dq1_0 = 0.05;
dq2_0 = 0;

Q1_MAX = pi+ pi/18;
Q2_MAX = pi/18;
diffmax = pi/12;


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

if use_lqr == 1
    
    [K,s,e] = lqr(A,b,Q,R);

    K_p =  K(1:2);
    K_d =  K(3:4);
else
    K_p = zeros(1,2);
    K_d = zeros(1,2);
end

%% exec


x = zeros(4,time);
dstate = zeros(4,time);

u = zeros(1,time);              %the final input
tau = zeros(1,time);            %LQR output

h = zeros(1,time);
q_des = zeros(2,time);
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
    
%     Q_DES(2) = Q2E + (pi/144)*sin(i*pi/125);
    
    q_des(:,i) = Q_DES;

    tau(i) = - K_p*(x(1:2,i)-Q_DES) - K_d*x(3:4,i);

    if use_cbf == 1
        [h(i), u(i)] =  CBFcontroller(x(:,i),tau(i)+ tau_eq, true);
%         fprintf("\n\n u: %f\n",(u(i)));
%         fprintf("tau: %f\n",(tau(i)));

    else
        [h(i),~] =  CBFcontroller(x(:,i),tau(i)+ tau_eq,false);
        u(i) = tau(i)+ tau_eq;
    end

    dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1],[x(:,i)',u(i)]);
    x(:,i+1) = x(:,i) + dt*dstate(:,i);

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
    N = (1 : time) * dt;

  figure;
        hold on
        plot(x(1,1),x(2,1),'o','Color','red');
        plot(x(1,1:time),x(2,1:time))
        plot(x(1,time),x(2,time),'o','Color','green')
       %ezplot(@(X,Y)(0.5*(diffmax^2 -(X- pi)^2 - Y^2)));
       th = 0:pi/50:2*pi;
        xunit = diffmax * cos(th) + pi;
        yunit = diffmax * sin(th) + 0;
       plot(xunit, yunit);
       xlim([pi-pi/10,pi+pi/10])
       ylim([-pi/10,pi/10])
    hold off

    
      figure
        plot(N,h)
        xlabel('time')
        ylabel('h')
        title('Control Barrier function')
        grid on
        
        
if show_plots
    figure();
    plot(N,x(2,1:time),N,q_des(2,1:time))
     yline(Q2_MAX,'--');
     yline(-Q2_MAX,'--');
    xlabel('s')
    ylabel('q2')
    legend('q2','q2_des');
    yline(pi/18)
    title('q2 angle')
    ylim([-0.5 0.4])

    figure
    plot(N,x(1,1:time),N,q_des(1,1:time))
%     yline(pi+Q1_MAX);
%     yline(pi-Q1_MAX);
    ylim([-3 3])
    xlabel('s')
    ylabel('q1')
    legend('q1','q1_des');
    yline(pi+pi/18)
    yline(pi-pi/18)
    title('q1 angle')

figure
plot(N,x(3,1:time),N,x(4,1:time));
    yline(1);
    yline(-1);
    xlabel('time');
    ylabel('dq');
    legend('dq1','dq2');
    title('dq velocity');
    ylim([-1.5 1.5])

   figure
    plot(x(2,1:time),x(4,1:time))
    yline(1);
    yline(-1);
    xlabel('q2')
    ylabel('dq2')
    title('q2 and dq2')
    ylim([-1.5 1.5])
    grid on

    if(use_cbf == 1)
        figure;
       figure
        plot(N,u(1,:),N,tau)
        xlabel('time')
        ylabel('input torque')
        legend('u','LQR')
        title('Final control input vs LQR component')

        figure
        plot(N,h)
        xlabel('time')
        ylabel('h')
        title('Control Barrier function')
        grid on
        
       
    else
        figure(5)
        plot(N,u(1,:))
        xlabel('time')
        ylabel('tau1')
        title('LQR input')
    end
end

%VideoSaver(x,time,l1,l2);
%VideoSaver genera un file .avi che deve essere accelarato 5x imo oppure
%ridurre il framerate


%% CBF
function [h, u] = CBFcontroller(x,u_des, active)
global a1 a2 a3 a4 a5 f1 f2 cbf grad f g alpha diffmax

syms q1 q2 dq1 dq2;

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

F = diag([f1,f2]);   %friction matrix

f = [dq1;
    dq2;
    -M\(c+e+F*[dq1;dq2])
    ];

g = [0,0;0,0; -M];
g = g(:,1);

alpha = 1;
par1 = 1;
par2 = 1;
%%%%
cbf = 0.5*(diffmax^2 - (q1-pi)^2 - q2^2 -par1*dq2^2);
grad = [-(q1-pi), -q2,  0,   -par1*dq2];

A = -grad*g;
b = grad*f + alpha*cbf;
H = 1;
% %%%%%
disp("A:"),
disp(A);
disp("B:"),
disp(b);

if ~active
    h = cbf;
    u = 0;
    return;
end
options = optimset('display','off');

h = cbf;
u = quadprog(H,-u_des,A,b,[],[],[],[],[],options);
end

% function [c,ceq] = mycon(u)
% global cbf grad f g alpha
% c = -grad*(f+g*u)-alpha*cbf;
% ceq =  []; 
% end