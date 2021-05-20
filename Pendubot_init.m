% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>

clc;
clear all;
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
I2 = 0.005;
G = 9.81;            % constante di accelerazione gravitazionale
f1 = 0.5;
f2 = 0.5;
F = diag([f1,f2]);   %friction matrix


time = 500;
dt = 0.005;
q1_0 = pi;
q2_0 = 0;
dq1_0 = 0;
dq2_0 = 0;



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

f =[dq1;
    dq2;
    -inv(M)*(n+F*[dq1;dq2])
    ];

g =[0,0;
    0,0;
    inv(M)];

dx=f+g*[u1;0];
%% LQR INIT

Q1E = pi;

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
     -inv(Me)*He,   -inv(Me)*F];

 disp("A: ");
 disp(A);
 
b = [0; 0; Me\n];

disp("b: ");
disp(b);

[K,s,e] = lqrd(A,b,Q,R,0.005);

K_p = -K(1:2);
K_d = -K(3:4);



%% exec
x = zeros(4,time);

dstate = zeros(4,time);
dstate(:,1) =subs(dx,[q1,q2,dq1,dq2,u1],[q1_0,q2_0,dq1_0,dq2_0,0]);
x(:,1) = [q1_0,q2_0,dq1_0,dq2_0]';

%%%%%%%%%%%%%%%%%%%%
fig = figure();
    j1x = l1*sin(x(1,1));
    j1y = -l1*cos(x(1,1));
    j2x = j1x + l2*sin(x(2,1)+x(1,1));
    j2y = j1y - l2*cos(x(2,1)+x(1,1));
    
    hold on
    joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
    joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
    tip = plot(j2x,j2y,'*','MarkerFaceColor','red','Color','red');
    xlim([-1 1]);
    ylim([-1,1]); 
    hold off 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disturb = 0

for i = 1:time
    tau = 0;
    tau = K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i)+ tau_eq;
    if i <3
        disturb = 10;
    else 
        disturb = 0;
    end
     
    u =  CBFcontroller(f,g,x(:,i),tau)+ [disturb;0];
   
    dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1],[x(:,i)',u(1)]);
    x(:,i+1) = x(:,i) + dt*dstate(:,i);
    %disp("x(:,i+1) = "+x(:,i)+" + "+dt+"*"+dstate(:,i)+" :");  
    %disp(x(:,i+1));
    t(i+1) = i*dt;
    figure(fig);
    delete(joint1);
        delete(joint2);
        delete(tip);
        j1x = l1*sin(x(1,i));
        j1y = -l1*cos(x(1,i));
        j2x = j1x + l2*sin(x(2,i)+x(1,i));
        j2y = j1y - l2*cos(x(2,i)+x(1,i));

        hold on
        joint1 = plot([0,j1x],[0,j1y],'-o','MarkerFaceColor','blue','Color','blue');
        joint2 = plot([j1x,j2x],[j1y,j2y],'-','MarkerFaceColor','blue','Color','blue');
        tip = plot(j2x,j2y,'*','MarkerFaceColor','red','Color','red');
        xlim([-1 1]);
        ylim([-1,1]);
        hold off
        drawnow();
end

close(fig);

%% plot
%PenduPlot(x,time,l1,l2);



function u = CBFcontroller(f,g,x, u_des)
% syms q1 q2 dq1 dq2 real
% f_attuale = vpa(subs(f,[q1,q2,dq1,dq2],x'));
% g_attuale = vpa(subs(g,[q1,q2,dq1,dq2],x'));
m1 = 20;              % giunto 1 mass
m2 = 20;              % giunto 2 mass
l1 = 0.5;             % lunghezza primo giunto
lc1 = 0.25;           % centro di massa giunto 1
l2 = 0.5;
lc2 = 0.25;
I1 = 0.05;
I2 = 0.05;
G = 9.81;            % constante di accelerazione gravitazionale
f1 = 0.5;
f2 = 0.5;

d1= l1-lc1;
d2 =l2-lc2;


M = [(l1^2 + 2*cos(x(2))*l1*lc2 + lc2^2)*m2 + m1*lc1^2 + I1 + I2, lc2*(lc2 + l1*cos(x(2)))*m2 + I2;
     (lc2^2 + l1*cos(x(2))*lc2)*m2 + I2,  m2*lc2^2 + I2];
 
c = [x(4)*l1*m2*sin(x(2))*(x(3) + x(4))*(d2 - l2) + x(3)*x(4)*l1*m2*sin(x(2))*(d2 - l2);
    -x(3)^2*l1*m2*sin(x(2))*(d2 - l2)];

e = [ G*(l1*m2*cos(x(1)) + lc1*m1*cos(x(1)) + lc2*m2*cos(x(1) + x(2)));
    G*lc2*m2*cos(x(1) + x(2))];
n = c+e;
F = diag([f1,f2]);   %friction matrix

f_attuale =[x(3);
    x(4);
    -inv(M)*(n+F*[x(3);x(4)])
    ];

g_attuale =[0,0;
    0,0;
    inv(M)];


disp("f_attuale: ");
disp(f_attuale);
disp("g_attuale: ");
disp(g_attuale);

%cbf
Q1E = pi;

q1e = pi;
q2e = pi - Q1E;

Q2E_MAX = pi/12;
c = 0.1;
alpha = 1;
cbf = 0.5*(Q2E_MAX^2-(x(2)-q2e)^2-c*x(4)^2);
grad = [q1e - x(1), 0, -c*x(3), 0];
disp("CBF: " );
disp(cbf);

Aeq = [];
beq = [];

H = eye(2);
f = -[u_des; 0];
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

u= quadprog(H,f,A,b,Aeq,beq,[],[],[u_des,0]',options);
end