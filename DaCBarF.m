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

h_dot_est = zeros(1,T);

q1_0 = pi+pi/12;
q2_0 = -pi/12;
dq1_0 = 0;
dq2_0 = 0;

x_0 = [q1_0, q2_0, dq1_0, dq2_0];

cbf = 0.5*(1-c_param*dq2^2);
h_dot_est_0 = 0.5*(1-2*c_param*dq2);

% w = ones(1,T); for now w = 1 for all j

for i = 1:T
    % Sample initial conditions
    if i == 0
        D_curr = experiment(x_0,h_dot_est_0,i); % Execute experiment
    else
        D_curr = experiment(x_0,h_dot_est(i),i);
    end
    D(:,time*(i-1)+1:time*i) = D_curr; %Aggregate dataset
    % Fit estimators
    % h_dot_est(i) = h_dot_est_0 + a*u + b % Update derivative estimators
    % Update controller will be directly in the experiment
end


function D = experiment(x_0,h_dot_est,num_exp)

    q1_0  = x_0(1);
    q2_0  = x_0(2);
    dq1_0 = x_0(3);
    dq2_0 = x_0(4);

    m2_var = 2;

    x = zeros(4,time);
    dstate = zeros(4,time);

    u = zeros(2,time);              %the final input
    tau = zeros(1,time);            %LQR output
    
    D = zeros(5,time);              %experiment output
    h = zeros(1,time);
    h_dot = zeros(1,time);
    q_des = zeros(2,time);
    dstate(:,1) = subs(dx,[q1,q2,dq1,dq2,u1,m2_unc],[q1_0,q2_0,dq1_0,dq2_0,0,m2_var]);
    x(:,1) = [q1_0,q2_0,dq1_0,dq2_0]';

    for i = 1:time
        q_des(:,i) = Q_DES;

        tau(i) = K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i) + tau_eq;
        
        if num_exp == 0 %first experiment executed with k0, CBF-QP controller
            [h(i), h_dot(i), u(1,i)] = CBFcontroller(x(:,i),tau(i));
        else
            [h(i), h_dot(i), u(1,i)] =  LCBFcontroller(x(:,i),tau(i),h_dot_est);
        end
        
        % Mass pertubation
        a = 90;
        b = 110;
        r = (b-a).*rand() + a;
        m2_var = 2*r/100;
        
        dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1,m2_unc],[x(:,i)',u(1,i),m2_var]);
        x(:,i+1) = x(:,i) + dt*dstate(:,i);
        
        D(1:4,i) = x(:,i);
        D(5,i) = h_dot(i);
        
    end
end

function [h, h_dot, u] = CBFcontroller(x,u_des)
    global a1 a2 a3 a4 a5 f1 f2
    syms q1 q2 dq1 dq2;

    q1  = x(1);
    q2  = x(2);
    dq1 = x(3);
    dq2 = x(4);

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

    cbf = 0.5*(1-c_param*dq2^2);
    cbf_dot = 0.5*(1-2*c_param*dq2);
    grad = [0,0,0,-c_param*dq2];

    H = 1;
    f_qp = -u_des;

    A = -grad*g;
    b = grad*f+alpha*cbf;

    % disp("A:"); disp(A);
    % disp("b:"); disp(b);
    options = optimset('display','off');
    h = cbf;
    h_dot = cbf_dot;
    u = quadprog(H,f_qp,A,b,[],[],[],[],[],options);
end

function [h,h_dot,u] = LCBFcontroller(x,u_des,h_dot_est)
    syms q1 q2 dq1 dq2;

    q1  = x(1);
    q2  = x(2);
    dq1 = x(3);
    dq2 = x(4);

    % Constraint depends on S_dot (derivative estimator)
    % S_dot = h_dot + a^Tu + b

    c_param = 1;
    alpha = 20;
    
    cbf = 0.5*(1-c_param*dq2^2);
    cbf_dot = 0.5*(1-2*c_param*dq2);
    grad = [0,0,0,-c_param*dq2];

    H = 1;
    f_qp = -u_des;

    % A = ...
    % b = ... + alpha*cbf;

    % disp("A:"); disp(A);
    % disp("b:"); disp(b);
    options = optimset('display','off');
    h = cbf;
    h_dot = cbf_dot;
    u = quadprog(H,f_qp,A,b,[],[],[],[],[],options);
end

% TODO
function erm()
end

% TODO
function augment()
end