% Cartpole implementation
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear all;
close;

%% configuration

use_cbf = true;
use_lqr = true;
disturbance = true;
show_animation = false;
show_plots = true;

%% param setting
syms x theta dx dtheta ddx ddtheta u1 real

mc = 5;              % giunto 1 mass
mp = 1;              % giunto 2 mass
l = 2;             % lunghezza  giunto
G = 9.81;            % constante di accelerazione gravitazionale

time = 1200;
dt = 0.01;
x_0 = 4.5;
theta_0 = pi-pi/12;
dx_0 = 0%3;
dtheta_0 =0 %-1;

THETAE_MAX = pi/4;%pi/12;

X_EQ = 0;
THETA_EQ = pi+pi/12;

xe = X_EQ;
THETA_LQR = pi;

Q_DES = [X_EQ; THETA_LQR];

%% model without uncertainties

M = [mc+mp,            mp*l*cos(theta);
     mp*l*cos(theta),    0.5+mp*l^2];
 

c = [0 -mp*l*dtheta*sin(theta);
     0              0;];

e = [0; mp*G*l*sin(theta)];

f =[dx;
    dtheta;
    -inv(M)*(c*[dx;dtheta]+e)
    ];

g =[0,0;
    0,0;
    inv(M)];

d_state=f+g*[u1;0];

%% LQR INIT


n = [1 0]';

Q = diag([1 1 10 100]);
R = 1;

tau_eq = mp*G*l*sin(THETA_LQR);

Me = [mc+mp,            mp*l*cos(THETA_LQR);
     mp*l*cos(THETA_LQR),    mp*l^2];

He = [0 0;
      0 mp*G*l;];

A = [zeros(2), eye(2);
     Me\He,  zeros(2)];
 
b = [0; 0; Me\n];

[K,s,e] = lqr(A,b,Q,R);

%% exec

state = zeros(4,time);
dstate = zeros(4,time);
u = zeros(1,time);
tau = zeros(1,time);
h = zeros(1,time);
q_des = zeros(2,time);
state(:,1) = [x_0,theta_0,dx_0,dtheta_0]';
dstate(:,1) =subs(d_state,[x,theta,dx,dtheta,u1],[state(:,1)',0]);

if(show_animation == 1)
    fig = figure();
    jx = l*sin(state(2,1));
    jy = -l*cos(state(2,1));

    hold on
    cart = plot(state(1,1),0,'gs','MarkerFaceColor','red','markers',22);
    joint= plot([state(1,1),state(1,1)+jx],[0,jy],'-','MarkerFaceColor','blue','Color','blue');
    pole = plot(state(1,1)+jx,jy,'-o','MarkerFaceColor','red','Color','red');
    xlim([state(1,1)-1 state(1,1)+1]);
    ylim([-0.5,2]);
    axis equal
    hold off 
end

%% Simulation loop
for i = 1:time
    Q_DES(2) = THETA_EQ;
    q_des(2,i) = Q_DES(2);
    tau(i) = -K*(state(:,i)-[Q_DES; 0; 0]) + tau_eq;
    
    if use_cbf == 1
        [h(i), u(i)] =  CBFcontroller(state(:,i),THETAE_MAX,THETA_LQR,tau(i));
%         disp(tau-u(1,i));
    else
        u(i) = tau(i);
    end
    u(i) = u(i);
    dstate(:,i) = subs(d_state,[x,theta,dx,dtheta,u1],[state(:,i)',u(i)]);
    state(:,i+1) = state(:,i) + dt*dstate(:,i);
    %disp("x(:,i+1) = "+x(:,i)+" + "+dt+"*"+dstate(:,i)+" :");  
    %disp(state(:,i+1));
    %t(i+1) = i*dt;

 if(show_animation == 1)
    figure(fig);
    jx = l*sin(state(2,i));
    jy = -l*cos(state(2,i));
    delete(cart);
    delete(joint);
    delete(pole);
    hold on
    cart = plot(state(1,i),0,'gs','MarkerFaceColor','red','markers',22);
    joint= plot([state(1,i),state(1,i)+jx],[0,jy],'-','MarkerFaceColor','blue','Color','blue');
    pole = plot(state(1,i)+jx,jy,'-o','MarkerFaceColor','red','Color','red');
    xlim([state(1,i)-1 state(1,i)+1]);
    ylim([-0.5,2]);
    axis equal
    hold off 
 end
end


%% plot
% show the state of the robot
if show_plots
    N = (1 : time) * dt;
    figure();
    %subplot(2,2,1);
    plot(N,state(2,1:time),N,q_des(2,1:time))
     yline(THETA_EQ+ THETAE_MAX,'--');
     yline(THETA_EQ-THETAE_MAX,'--');
    xlabel('time')
    ylabel('theta')
    legend('theta','theta_des');
    title('theta angle')
    ylim([-2.7 3.8])


    %subplot(2,2,2);
    figure;
    plot(N,state(1,1:time),N,q_des(1,1:time))
%     yline(pi+Q1_MAX);
%     yline(pi-Q1_MAX);
    ylim([-3 3])
    xlabel('time')
    ylabel('x')
    legend('x','x_des');
    title('x dynamics')

   % subplot(2,2,3)
%    figure
%     plot(N,state(3,1:time),N,state(4,1:time));
%     yline(1);
%     yline(-1);
%     xlabel('dtheta');
%     ylabel('dx');
%     legend('dx','dtheta');
%     title('dq velocity');
%     ylim([-1.5 1.5])

    %subplot(2,2,4)
    figure
    hold on
    plot(state(2,1:time),state(4,1:time))
    ezplot(@(X,Y)(0.5*(THETAE_MAX^2-(X-THETA_LQR)^2-Y^2)))
    hold off
    xlabel('theta')
    ylabel('dtheta')
   % xline(THETA_EQ+ THETAE_MAX,'--');
    %xline(THETA_EQ-THETAE_MAX,'--');
    title('theta and dtheta')
    ylim([-1.2 1.2])
    xlim([2.5 3.7])
   
    grid on

    if(use_cbf == 1)
        figure;
        %subplot(1,2,1)
        plot(N,u(1,:),N,tau)
        xlabel('time')
        ylabel('input torque')
        legend('u','LQR')
        title('Final control input vs LQR component')

       % subplot(1,2,2)
       figure;
        plot(N,h)
        xlabel('time')
        ylabel('h')
        title('Control Barrier function')
        grid on
    else
        figure(5)
        plot(N,u(1,:),N,tau)
        xlabel('time')
        ylabel('input torque')
        legend('u','LQR')
        title('Final control input vs LQR component')
    end
end

function [h, u] = CBFcontroller(state,THETAE_MAX,thetae,u_des)

mc = 10;              % giunto 1 mass
mp = 1;              % giunto 2 mass
l = 1;             % lunghezza primo giunto
G = 9.81;            % constante di accelerazione gravitazionale

theta = state(2);
dx = state(3);
dtheta = state(4);

M = [mc+mp,            mp*l*cos(theta);
     mp*l*cos(theta),    mp*l^2]; 

c = [0 -mp*l*dtheta*sin(theta);
     0              0;];

e = [0; mp*G*l*sin(theta)];

n = [1 0]';

f =[dx;
    dtheta;
    -M\(c*[dx;dtheta]+e)
    ];

g =[0;
    0;
    M\n];

%cbf
c = 1;
cbf = 0.5*(THETAE_MAX^2-(theta-thetae)^2-c*dtheta^2);
grad = [0, thetae - theta, 0, -c*dtheta];

H = 1;
f_qp = -u_des;
A = -grad*g;
b = grad*f+cbf;
options = optimset('display','off');

h = cbf;
u = quadprog(H,f_qp,A,b,[],[],[],[],[],options);
end