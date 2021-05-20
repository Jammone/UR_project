% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear;
close;
%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 u1 real

m1 = 20;              % giunto 1 mass
m2 = 20;              % giunto 2 mass
l1 = 0.5;             % lunghezza primo giunto
lc1 = 0.25;           % centro di massa giunto 1
l2 = 0.5;
lc2 = 0.25;
I1 = 0.05;
I2 = 0.05;

g = 9.81;            % constante di accelerazione gravitazionale

PARAM= [m1,m2,I1,lc1,l2,lc2,I1,I2];
INITIAL_CONDITIONS = [-pi/2,0,0,0];
KP_GAIN = [10,0;
            0,10];
KD_GAIN = [9,0;
            0,9];
DES_POS = [pi/2,0];
%% model without uncertainties
d1= l1-lc1;
d2 =l2-lc2;
% M*[ddq1; ddq2] + n = [u;0]

% a1 = m2*l1^2+m2*lc2^2 + I2 +I1+ m1*lc1^2
% a2= m2*l1*lc2;
% a3 =m2*lc2^2+I2;

    syms m1 m2 I1 I2 lc1 lc2 l1 l2 real
   
M = [(l1^2 + 2*cos(q2)*l1*lc2 + lc2^2)*m2 + m1*lc1^2 + I1 + I2, lc2*(lc2 + l1*cos(q2))*m2 + I2;
     (lc2^2 + l1*cos(q2)*lc2)*m2 + I2,  m2*lc2^2 + I2];
 
c = [dq2*l1*m2*sin(q2)*(dq1 + dq2)*(d2 - l2) + dq1*dq2*l1*m2*sin(q2)*(d2 - l2);
    -dq1^2*l1*m2*sin(q2)*(d2 - l2)];

e = [ g*(l1*m2*cos(q1) + lc1*m1*cos(q1) + lc2*m2*cos(q1 + q2));
    g*lc2*m2*cos(q1 + q2)];
n = c+e;

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
% %plot gragico dinamico avanzato
 fig = figure();
 j1x = l1*cos(INITIAL_CONDITIONS(1));
 j1y = l1*sin(INITIAL_CONDITIONS(1));
 j2x = j1x + l2*cos(INITIAL_CONDITIONS(2)+INITIAL_CONDITIONS(1));
 j2y = j1y + l2*sin(INITIAL_CONDITIONS(2)+INITIAL_CONDITIONS(1));
 
 hold on
 joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
 joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
 tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
 xlim([-1 1]);
 ylim([-1,1]);
 p = plot(0,0,'s','MarkerFaceColor','blue' ); 
%  
hold off 
 
 for k = 2:length(Sq1.Data)
     figure(fig)
     delete(joint1);
     delete(joint2);
     delete(tip);
     j1x = l1*cos(Sq1.Data(k));
     j1y = l1*sin(Sq1.Data(k));
     j2x = j1x + l2*cos(Sq2.Data(k)+Sq1.Data(k));
     j2y = j1y + l2*sin(Sq2.Data(k)+Sq1.Data(k));
 
 hold on
 joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
 joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
 tip = plot(j2x,j2y,'>','MarkerFaceColor','red','Color','red');
 xlim([-1 1]);
 ylim([-1,1]);
 hold off
     pause(0.005);%faster than drawnow
     
 end





%% Optimal control
 %J= 0.5*norm(u-Kp(-state))^2;
 
 


