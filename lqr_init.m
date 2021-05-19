clc;
clear;
close;

syms q1 q1dot q2 q2dot;
x = [q1 q1dot q2 q2dot]';

m1 = 1;            % giunto 1 mass
m2 = 1;              % giunto 2 mass
l1 = 0.5;             % lunghezza primo giunto
lc1 = l1/2;           % centro di massa giunto 1
l2 = 0.5;
lc2 = l2/2;
I1 = 0.1;
I2 = 0.1;

f1 = 0.5;
f2 = 0.5;

F = diag([f1,f2]);

g = 9.81;            % constante di accelerazione gravitazionale

% PARAM = [m1,m2,I1,lc1,l2,lc2,I1,I2];

a1 = I1+m1*lc1^2+I2+m2*(l1^2+lc2^2);
a2 = m2*l1*lc2;
a3 = I2+m2*lc2^2;
a4 = g*(m1*lc1+m2*l1);
a5 = g*m2*lc2;

COEFF = [a1 a2 a3 a4 a5 f1 f2];

INITIAL_CONDITIONS = [pi + pi/12, -pi/12, 0, 0];
Q1E = pi-pi/12;

q1e = Q1E;
q2e = pi - Q1E;

Q_DES = [q1e; q2e];

n = [1 0]';

Q = diag([10 10 1 1]);
R = 1;

tau_eq = [a4*sin(q1e)+a5*sin(q1e+q2e)];

Me = [a1 + 2*a2*cos(q2e), a3+a2*cos(q2e);
         a3+a2*cos(q2e),          a3];

He = [a4*cos(q1e)+a5*cos(q1e+q2e),  a5*cos(q1e+q2e);
           a5*cos(q1e+q2e),         a5*cos(q1e+q2e)];
       
A = [zeros(2), eye(2);
     -Me\He,   -Me\F];

b = [0; 0; Me\n];

[K,s,e] = lqrd(A,b,Q,R,0.005);

K_p = -K(1:2);
K_d = -K(3:4);

% CBF definition
Q1E_MAX = pi/12;
c = 0.1;
h = 0.5*(Q1E_MAX^2-(q1-q1e)-c*q1dot);

% QP definition
lb = -100;
ub = 100;

sim("lqr_pendubot.slx");
%plot gragico dinamico avanzato
 fig = figure();
 j1x = l1*sin(INITIAL_CONDITIONS(1));
 j1y = -l1*cos(INITIAL_CONDITIONS(1));
 j2x = j1x + l2*sin(INITIAL_CONDITIONS(2)+INITIAL_CONDITIONS(1));
 j2y = j1y - l2*cos(INITIAL_CONDITIONS(2)+INITIAL_CONDITIONS(1));
 
 hold on
 joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
 joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
 tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
 xlim([-1 1]);
 ylim([-1,1]);
 p = plot(0,0,'s','MarkerFaceColor','blue' ); 
 
hold off 
 
 for k = 2:length(Sq1.Data)
     figure(fig)
     delete(joint1);
     delete(joint2);
     delete(tip);
     j1x = l1*sin(Sq1.Data(k));
     j1y = -l1*cos(Sq1.Data(k));
     j2x = j1x + l2*sin(Sq2.Data(k)+Sq1.Data(k));
     j2y = j1y - l2*cos(Sq2.Data(k)+Sq1.Data(k));
 
 hold on
 joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
 joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
 tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
 xlim([-1 1]);
 ylim([-1,1]);
 hold off
     pause(0.005);%faster than drawnow
     
 end
