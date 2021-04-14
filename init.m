% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>


%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 u1 real

m1 = 5;              % giunto 1 mass
m2 = 5;              % giunto 2 mass
l1 = 0.5;             % lunghezza primo giunto
lc1 = 0.25;           % centro di massa giunto 1
l2 = 0.5;
lc2 = 0.25;
I1 = 0.05;
I2 = 0.05;

g = 9.81;            % constante di accelerazione gravitazionale

PARAM= [m1,m2,I1,lc1,l2,lc2,I1,I2];

%% model without uncertainties
% M*[ddq1; ddq2] + n = [u;0]
a1 = I1+ I2+ m1*lc1^2+m2*(l1^2+lc2^2);
a2 = m2*l1*lc2^2;
a3 = I2 + m2*lc2^2;
a4 = g*(m1*lc1+m2*l1);
a5 = g*m2*lc2;

M = [a1+2*a2*cos(q2), a3+a2*cos(q2);
     a3+a2*cos(q2), a3];
 
c = [-a2*dq2*(dq2+2*dq1)*sin(q2);
    a2*dq1^2*sin(q2)];
e = [ a4*sin(q1) + a5*sin(q1+q2);
    a5*sin(q1+q2)];
n = c+e;
% dstate = f(state) + g(state)*u
state = [q1,q2,dq1,dq2]';

f = [dq1;
    dq2;
    -inv(M)*n(1);
    -inv(M)*n(2)];

g =[0,0;
    0,0;
    inv(M);
    0,0];

dstate = f + g*[u1;0];
%% Optimal control
 J= 0.5*norm(u-Kp(-state))^2;
 
 

