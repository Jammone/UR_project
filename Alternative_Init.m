clc;
clear;
close;
%% param setting
syms q1 q2 dq1 dq2 ddq1 ddq2 m1 m2 I1 I2 l1 l2 d1 d2 lc1 lc2 g real



pc1 = [-d1;0;0];
pc2 = [-d2;0;0];
q=[q1,q2]; dq=[dq1,dq2];



DH_TABLE = [      0,        l1,      0,        q1;
                  0,        l2,     0,         q2;
            ];

                
%% Link 1 setup

I11 = [0,0,I1];
I22 = [0,0,I2];


L1   = Link('revolute','a',DH_TABLE(1,2));

L1.m = m1;                               %mass 
L1.I = eye(3)*I11';         %Inertia matrix, diagonal here
L1.r = pc1;                                                                  

%% Link 2 setup
 L2   = Link('revolute','a',DH_TABLE(2,2));
 L2.m = m2;
 L2.I = eye(3)*I22'; 
 L2.r = pc2;
 
 %% Creation of the robot
bot = SerialLink([L1,L2],'gravity',[0,+g,0]);     %feed him with the links

M = vpa(bot.inertia(q))       %inertia matrix, sometimes weid very little numbers
% might appear, just rewrite it without the 0.123456789e-33*element stuff
S = bot.coriolis(q,dq)
C = S*dq'
G= bot.gravload(q)'

PARAM= [m1,m2,I1,lc1,l2,lc2,I1,I2];

NPARAM = [M,C,G];

INITIAL_CONDITIONS = [0,0,0,0];
KP_GAIN = [10,0;0,10];
KD_GAIN = [2,0;0,9];
DES_POS = [pi/2,0];
M= subs(M,[d1,d2],[l1-lc1,l2-lc2]);
M= simplify(M)
G= subs(G,[d1,d2],[l1-lc1,l2-lc2]);
G=simplify(G)
%sim("Modello_Simulink.slx");
