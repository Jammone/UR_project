% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear all;
close;

%% param setting
syms x phi dx dphi ddx ddphi u real

time = 10;
dt = 0.01;
time = time/dt;
x_0 = 0;
phi_0 = pi/18;
dx_0 = 0;
dphi_0 = 0;


%% model without uncertainties

f =[dx;
    dphi;
    (cos(phi)*(11.5*dx+9.8*sin(phi))+68.4*dx-1.2*dphi^2*sin(phi))/(cos(phi)-24.7);
    (-58.8*dx*cos(phi)-243.5*dx-sin(phi)*(208.3+dphi^2*cos(phi)))/(cos(phi)^2-24.7);];

g =[0;
    0;
    (-1.8*cos(phi)-10.9)/(cos(phi)-24.7);
    (9.3*cos(phi)+38.6)/(cos(phi)^2-24.7)];

d_state=f+g*u;

%% PD control gain
K_p = 1;
K_d = 0.1;

phi_des = 0;

%% exec
state = zeros(4,time);
dstate = zeros(4,time);
input = zeros(1,time);
dstate(:,1) = subs(d_state,[x,phi,dx,dphi,u],[x_0,phi_0,dx_0,dphi_0,0]);
state(:,1) = [x_0,phi_0,dx_0,dphi_0]';

%% Show initial state
l1  = 1;
fig = figure();
j1x = -l1*sin(state(2,1));
j1y = l1*cos(state(2,1));

hold on
wheel = plot(state(1,1),0,'O','MarkerFaceColor','red','Color','red');
joint1 = plot([state(1,1),j1x],[0,j1y],'-','MarkerFaceColor','blue','Color','blue');
axis equal

xlim([-2 2]);
ylim([-2,2]); 
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for i = 1:time
    input(1,i) = K_p*(phi_des - state(2,i)) - K_d*state(4,i);
%     disp(input(1,i));
    dstate(:,i) = subs(d_state,[x,phi,dx,dphi,u],[state(:,i)', input(1,i)]);
    state(:,i+1) = state(:,i) + dt*dstate(:,i);

%     disp("x(:,i+1) = "+state(:,i)+" + "+dt+"*"+dstate(:,i)+" :");
%     disp(x(:,i+1));
%     t(i+1) = i*dt;
    
    fig = figure(fig);
    delete(wheel);
    delete(joint1);
    j1x = -l1*sin(state(2,i+1));
    j1y = l1*cos(state(2,i+1));
    disp(j1x^2+j1y^2);

    hold on
    wheel = plot(state(1,i+1),0,'O','MarkerFaceColor','red','Color','red');
    joint1 = plot([state(1,i+1),j1x],[0,j1y],'-','MarkerFaceColor','blue','Color','blue');
    xlim([-2,2]);
    ylim([-2,2]); 
%   hold off 
    drawnow();
end

close(fig);

%% plot
%PenduPlot(x,time,l1,l2);

% show the state of the robot

N = (1 : i-1) * dt;

figure(1)
plot(N,state(1,1:i-1),N,state(2,1:i-1));
xlabel('time');
ylabel('')
legend('x','phi');
title('x & phi')

figure(2)
plot(N,state(3,1:i-1),N,state(4,1:i-1));
xlabel('time');
ylabel('');
legend('dx','dphi');
title('velocity');

figure(3)
plot(state(1,1:i-1),state(3,1:i-1));
xlabel('x');
ylabel('dx')
title('x & dx')

figure(4)
plot(state(2,1:i-1),state(4,1:i-1));
xlabel('phi');
ylabel('dphi')
title('phi & dphi')

figure(5)
plot(N,input(1,1:i-1));
xlabel('time');
ylabel('tau1');
title('input torque');