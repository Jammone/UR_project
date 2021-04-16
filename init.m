% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear;
close;
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
INITIAL_CONDITIONS = [pi,0,0,0];

Q = diag([1,1,10,10]);
R = [10];
EQ_POS = [pi,0];

q1e = EQ_POS(1);
q2e = EQ_POS(2);

s = [1 0]';

%% model without uncertainties
% M*[ddq1; ddq2] + n = [u;0]
a1 = I1+ I2+ m1*lc1^2+m2*(l1^2+lc2^2);
a2 = m2*l1*lc2;
a3 = I2 + m2*lc2^2;
a4 = g*(m1*lc1+m2*l1);
a5 = g*m2*lc2;

M = [a1+2*a2*cos(q2), a3+a2*cos(q2);
     a3+a2*cos(q2), a3];
 
c = [a2*dq2*(dq2+2*dq1)*sin(q2);
    a2*dq1^2*sin(q2)];
e = [ a4*sin(q1) + a5*sin(q1+q2);
    a5*sin(q1+q2)];
n = c+e;

%Calculate He and Me and tau_eq
He = [a4*cos(q1e)+a5*cos(q1e+q2e), a5*cos(q1e+q2e);
      a5*cos(q1e+q2e),             a5*cos(q1e+q2e);];
Me = [a1+2*a2*cos(q2e), a3+a2*cos(q2e);
     a3+a2*cos(q2e), a3];
tau_eq =  a4*sin(q1e) + a5*sin(q1e+q2e);
%tau_eq = 0;

%Define approximation of A and B in function of He and Me
A = [zeros(2), eye(2);
     -inv(Me)*He, zeros(2)];
B = [0; 0; inv(Me)*s];
K = lqr(A,B,Q,R); %--> K will be a 1 x 4 vector
K_p = K(:,1:2);
K_d = K(:,3:4);

% dstate = f(state) + g(state)*u
state = [q1,q2,dq1,dq2]';

f = [dq1;
    dq2;
    -inv(M)*n];

g =[0,0;
    0,0;
    inv(M)];

dstate = f + g*[u1;0];

sim("Modello_Simulink.slx");
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
% 
% % Optimal control
%  J= 0.5*norm(u-Kp(-state))^2;