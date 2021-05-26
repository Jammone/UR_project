%% Segway

clc;
clear all;
close all;

syms mc mp l g A B Q R xd x theta d_x d_theta u1

mc = 10;  % Mass of the cart (wheel)
mp = 1;   % Mass of the pendulum (steer)
l = 1;    % Lenght of the link
g0 = 9.81;

% T0 = 0;
dt = 0.005; %1e-1;
time = 500;
% T = 20;
% tspan = T0:dt:T;


x0 = [-2 pi+1.1 0 0]';
xd = [0 pi 0 0 ]';

M = [mc+ mp -mp*l;
    -mp*l mp*l^2];

b = [1;0];

dtdq = [0 0;0 mp*g0*l];

A = [ zeros(2) eye(2);
    M\dtdq zeros(2)];

B = [zeros(2,1);M\b];

Q = diag([10 10 1 1]);
R = 1;


%% model without uncertainties

p = mc+mp*sin(theta)^2;
f = [d_x;
    d_theta;
    1/p*mp*sin(theta)*(l*d_theta^2+g0*cos(theta));
    1/(l*p)*(mp*l*d_theta^2*cos(theta)*sin(theta)-(mc+mp)*g0*sin(theta))];

g = [0;
    0;
    1/p;
    -1/(p*l)*cos(theta)];


dyn = f+g*u1;

%% LQR

% compute control u
x_e=x0; %-xd;
[K,s,e] = lqrd(A,B,Q,R,dt); %[K] = lqr(A,B,Q,R);
u = -K*x_e;

% compute dynamics
q = x0(2); dx = x0(3);dq = x0(4);
Me = [mc + mp  mp*l*cos(q); mp*l*cos(q) mp*l^2];

C = [0 -mp*l*dq*sin(q); 0 0];
G = [0 ;mp*g0*l*sin(q)];
dq = [dx;dq;];
B1 = [1;0];

ddq = Me\(B1*u-C*dq-G);
xdot = [dq;ddq];

%% exec
x = zeros(4,time);

dstate = zeros(4,time);
%dstate(:,1) = xdot;
dstate(:,1) = subs(dyn,[theta, d_x, d_theta,u1],[x0(2),x0(3),x0(4),u]);
x(:,1) = x0;

%% %%%%%%%%%%%%%%%%%%%%
fig = figure();
axp = l*cos(x(2,1))+x(1,1);
ayp = l*sin(x(2,1));

hold on
baseplot = plot(x(1,1),0,'gs','MarkerFaceColor','red','markers',22);
axeplot = plot([x(1,1),axp],[0,ayp],'-','MarkerFaceColor','blue','Color','blue');
xlim([-5 5]);
ylim([-5,5]);  
hold off 
% %%%%%%%%%%%%%%%%%%%%

for i = 1:time
    %u_ = -K*(x(:,i)-xd);
    u_ = -K*x(:,i);
    u =  CBFcontroller(f,g,x(:,i),u_);

    dstate(:,i) = subs(dyn,[theta, d_x, d_theta,u1],[x(2:4,i)',u]);
    x(:,i+1) = x(:,i) + dt*dstate(:,i);
    
    figure(fig);
    delete(baseplot);
    delete(axeplot);
   
    axp = l*cos(sin(x(2,i)))+x(1,i);
    ayp = l*sin(sin(x(2,i)));

    hold on
    baseplot = plot(x(1,i),0,'gs','MarkerFaceColor','red','markers',22);
    axeplot = plot([x(1,i),axp],[0,ayp],'-','MarkerFaceColor','blue','Color','blue');
    xlim([-5 5]);
    ylim([-5,5]);
    hold off
    drawnow();     
end
close(fig);


function u = CBFcontroller(f,g,x, u_des)

mc = 10;  % Mass of the cart (wheel)
mp = 1;   % Mass of the pendulum (steer)
l = 1;    % Lenght of the link
g0 = 9.81;

p = mc+mp*sin(x(2))^2;
f_attuale = [x(3);
             x(4);
             1/p*mp*sin(x(2))*(l*x(4)^2+g0*cos(x(2)));
             1/(l*p)*(mp*l*x(4)^2*cos(x(2))*sin(x(2))-(mc+mp)*g0*sin(x(2)))];

g_attuale = [0;
             0;
             1/p;
             -1/(p*l)*cos(x(2))];

disp("f_attuale: ");
disp(f_attuale);
disp("g_attuale: ");
disp(g_attuale);

% xd = [0 pi 0 0 ]';

xe = 0;
thetae = pi;
theta_max = pi-0.2;

c = 0.1;
alpha = 1;
cbf = 0.5*(theta_max^2-(x(2)-thetae)^2-c*x(4)^2);

grad = [0, thetae - x(2), 0, -c*x(4)];
disp("CBF: " );
disp(cbf);

H = 1.0; %eye(2);
f = -u_des; % -[u_des; 0];
A = -grad*g_attuale;
disp("grad ");
disp(grad);
disp("g_attuale " );
disp(g_attuale);

disp("A: " );
disp(A);
b = grad*f_attuale+alpha*cbf;
disp(b)
options = optimoptions('quadprog','Algorithm','active-set');

u= quadprog(H,f,A,b,[],[],[],[],u_des,options);
end