clc;
clear all; close all;

dt = 0.02;
sim_t = 5;

params.m1 = 10;
params.m2 = 10;
params.L1 = 1;
params.Lc1 = 0.5;
params.Lc2 = 0.5;
params.L2 = 1;
params.grav = 9.81;
params.I1 = 0.005;
params.I2 = 0.005;
params.cbf.rate = 1;
params.weight= [1,1];

pendulum = Pendubot(params);

odeFun = @pendulum.dynamics;
controller = @pendulum.ctrlCbfQp;
odeSolver = @ode45;

x0 =[pi/2, 3,0,0]';

%% TOFIX
total_k = ceil(sim_t / dt);
x = x0;
t = 0;   
% initialize traces.
xs = zeros(total_k, pendulum.xdim);

ts = zeros(total_k, 1);

us = zeros(total_k-1, 2);

Vs = zeros(total_k-1, 1);
xs(1, :) = x0';
ts(1) = t;
disp("hi");

%% inizializzo il grafico
 %plot gragico dinamico avanzato
 fig = figure();
 j1x = pendulum.L1*cos(x0(1));
 j1y = pendulum.L1*sin(x0(1));
 j2x = j1x + pendulum.L2*cos(x0(2)+x0(1));
 j2y = j1y + pendulum.L2*sin(x0(2)+x0(1));
 
 hold on
 joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
 joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
 tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
 xlim([-2.3 2.3]);
 ylim([-2.3,2.3]);
 p = plot(0,0,'s','MarkerFaceColor','blue' );  
hold off 
 
 
 
%% calcolo e simulo
for k = 1:total_k-1
    % Determine control input.
    % dV_hat: analytic Vdot based on model.
    %[u, B, feas, comp_time] = controller(x);        
    u=[0,0]';
    us(k, :) = u';
    % Run one time step propagation.
    [ts_temp, xs_temp] = odeSolver(@(t, s) odeFun(t, s, u), [t t+dt], x);
    
    x = xs_temp(end, :)';

    xs(k+1, :) = x';
    ts(k+1) = ts_temp(end);
    t = t + dt;
    
    
    figure(fig)
     delete(joint1);
     delete(joint2);
     delete(tip);
     j1x = pendulum.L1*cos(x(1));
     j1y = pendulum.L1*sin(x(1));
     j2x = j1x + pendulum.L2*cos(x(2)+x(1));
     j2y = j1y + pendulum.L2*sin(x(2)+x(1));
 
     hold on
     joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
     joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
     tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
    xlim([-2.3 2.3]);
    ylim([-2.3,2.3]);
     hold off
    
end

 


