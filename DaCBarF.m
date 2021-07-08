% DaCBarF
clc;
clear all;
close;

%% Pendubot settings
syms q1 q2 dq1 dq2 ddq1 ddq2 u1 m2_unc real
global Q_DES Q1E Q2E
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

%% LQR INIT
   
Q1E = pi-pi/18;

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

[K,s,e] = lqr(A,b,Q,R);

K_p_0 = -K(1:2);
K_d_0 = -K(3:4);

%% Algorithm settings

T = 10;
time = 1000;

D = zeros(4,T*time); % Initialize data set

q1_0 = pi+pi/12;
q2_0 = -pi/12;
dq1_0 = 0;
dq2_0 = 0;

x_0 = [q1_0, q2_0, dq1_0, dq2_0];

cbf = 0.5*(1-c_param*dq2^2);

grad = [0,0,0,-c_param*dq2]; %???
h_0 = grad*f + grad*g;

w = ones(1,T);

for i = 1:T
    % Sample initial conditions
    % D_curr = experiment % Execute experiment 
    % TODO: Implementation of simulator
    D(:,time*(i-1)+1:time*i) = D_curr; %Aggregate dataset
    % Fit estimators
    % Update derivative estimators
    % Update controller
end


function x = experiment(x_0,k)

    q1_0  = x_0(1);
    q2_0  = x_0(2);
    dq1_0 = x_0(3);
    dq2_0 = x_0(4);

    m2_var = 2;

    x = zeros(4,time);
    dstate = zeros(4,time);

    u = zeros(2,time);              %the final input
    tau = zeros(1,time);            %LQR output

    h = zeros(1,time);
    q_des = zeros(2,time);
    dstate(:,1) = subs(dx,[q1,q2,dq1,dq2,u1,m2_unc],[q1_0,q2_0,dq1_0,dq2_0,0,m2_var]);
    x(:,1) = [q1_0,q2_0,dq1_0,dq2_0]';

    for i = 1:time
        q_des(:,i) = Q_DES;

        tau(i) = K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i) + tau_eq;
        
        [h(i), u(1,i)] =  LCBFcontroller(x(:,i),tau(i));
        
        a = 90;
        b = 110;
        r = (b-a).*rand() + a;
        m2_var = 2*r/100;
        
        dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1,m2_unc],[x(:,i)',u(1,i),m2_var]);
        x(:,i+1) = x(:,i) + dt*dstate(:,i);
        
    end
end

function [h,u] = LCBFController(x,u_des)
    global a1 a2 a3 a4 a5 f1 f2 Q1E Q2E Q1_MAX Q2_MAX
    syms q1 q2 dq1 dq2;

    q1  = x(1);
    q2  = x(2);
    dq1 = x(3);
    dq2 = x(4);

%     n = [1 0]';
% 
%     M = [a1 + 2*a2*cos(q2), a3+a2*cos(q2);
%              a3+a2*cos(q2),          a3];
% 
%     c = [a2*sin(q2)*dq2*(dq2+2*dq1);
%              a2*sin(q2)*dq1^2];
% 
%     e = [a4*sin(q1)+a5*sin(q1+q2);
%                 a5*sin(q1+q2)];
% 
%     F = diag([f1,f2]);   %friction matrix
% 
%     f = [dq1;
%         dq2;
%         -M\(c+e+F*[dq1;dq2])
%         ];
% 
%     g = [0; 0; M\n];

% Constraint depends on S_dot (derivative estimator)
% S_dot = h_dot + a^Tu + b

    c_param = 1;
    alpha = 20;
    
    cbf = 0.5*(1-c_param*dq2^2);
    grad = [0,0,0,-c_param*dq2];

    H = 1;
    f_qp = -u_des;

    % A = ...
    % b = ... + alpha*cbf;

    % disp("A:"); disp(A);
    % disp("b:"); disp(b);
    options = optimset('display','off');
    h = cbf;
    u = quadprog(H,f_qp,A,b,[],[],[],[],[],options);
end

% TODO
function erm()
end

% TODO
function augment()
end