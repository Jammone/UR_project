% inverted pendulum (segway)
% progetto di <studente 1>, <studente 2> e <studente 3>


%% param setting
syms theta dtheta real

mb = 5;              % bar mass
MB = 30;             % body mass
l = 0.5;             % lunghezza asta 
g = 9.81;            % constante di accelerazione gravitazionale
I = 0.05;            % inerzia 


%% model without uncertainties
% M*[ddq1; ddq2] + n = [tau_1;0]

M = [mb+Mb, m*l*cos(theta);
     m*l*cos(theta), I+ mb*l*2;
     ];
n = [-mb*l*sin(theta)*dtheta;
    -mb+g*l*sin(theta)];

