% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>


%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 real

mb = 5;              % bar mass
MB = 30;             % body mass
l = 0.5;             % lunghezza asta 
g = 9.81;            % constante di accelerazione gravitazionale
I = 0.05;            % inerzia 


%% model without uncertainties
% M*[ddq1; ddq2] + n = [u;0]

M = [mb+Mb, m*l*cos(q2);
     m*l*cos(theta), I+ mb*l^2;
     ];
n = [-mb*l*sin(q2)*dq2;
    -mb+g*l*sin(q2)];

% dstate = f(state) + g(state)*u
state = [q1,q2,dq1,dq2]';

f = [dq1;
    dq2;
    inv(M)*(mb*l*sin(q2)*dq2);
    inv(M)*(mb*l*sin(q2))];

g =[0;
    0;
    inv(M);
    0];

dstate = f + g*[u;0];

