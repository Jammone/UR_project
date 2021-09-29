% DaCBarF
clc;
clear all;
close;

%% Pendubot settings

global Q_DES Q1E Q2E
global a1 a2 a3 a4 a5 f1 f2 K_p K_d tau_eq
global l1 lc1 l2 lc2 I1 I2 G F cbf_values exp_num
global time show_animation show_plots h_prev

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

F = diag([f1,f2]);   %friction matrix

%% model without uncertainties
a1 = I1+m1*lc1^2+I2+m2*(l1^2+lc2^2);
a2 = m2*l1*lc2;
a3 = I2+m2*lc2^2;
a4 = G*(m1*lc1+m2*l1);
a5 = G*m2*lc2;

%% LQR INIT

Q1E = pi+pi/18;

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

K_p = -K(1:2);
K_d = -K(3:4);

%% Algorithm settings

T = 7;
time = 1000;
cbf_values = zeros(T,time);
exp_num = 1;
show_animation = false;
show_plots = true;

h_dot_est = zeros(1,T);

q1_0 = pi-pi/12;
q2_0 = pi/12;
dq1_0 = 0;
dq2_0 = 0;

x_0 = [q1_0, q2_0, dq1_0, dq2_0];

a_est_net = [];
b_est_net = [];

% w = ones(1,T); for now w = 1 for all j

%% Algorithm execution

for i = 1:T
    % Sample initial conditions
    h_prev = 0.5;
    exp_num =  i;
    disp([num2str(i),'th experiment started'])
    D_curr = experiment(x_0,i,a_est_net,b_est_net); % Execute experiment
    disp('Learning started')
    D(:,time*(i-1)+1:time*i) = D_curr; %Aggregate dataset
    % Fit estimators
    [a_est_net, b_est_net] = erm(D);
    % Update controller will be directly in the experiment
    disp([num2str(i),'th experiment done'])
end

% Aux functions

function D = experiment(x_0,num_exp,a_est_net,b_est_net)

    syms q1 q2 dq1 dq2 ddq1 ddq2 u1 m1_unc m2_unc real
    global Q_DES l1 lc1 l2 lc2 I1 I2 G F K_p K_d tau_eq 
    global time dt show_animation show_plots cbf_values exp_num

    dt = 0.005;
    
    % Actual model with mass uncertainties
    
    a1_unc = I1+m1_unc*lc1^2+I2+m2_unc*(l1^2+lc2^2);
    a2_unc = m2_unc*l1*lc2;
    a3_unc = I2+m2_unc*lc2^2;
    a4_unc = G*(m1_unc*lc1+m2_unc*l1);
    a5_unc = G*m2_unc*lc2;

    M_unc = [a1_unc + 2*a2_unc*cos(q2), a3_unc+a2_unc*cos(q2);
             a3_unc+a2_unc*cos(q2),          a3_unc];


    c_unc = [a2_unc*sin(q2)*dq2*(dq2+2*dq1);
             a2_unc*sin(q2)*dq1^2];

    e_unc = [a4_unc*sin(q1)+a5_unc*sin(q1+q2);
                a5_unc*sin(q1+q2)];

    n_unc = c_unc+e_unc;


    f =[dq1;
        dq2;
        -inv(M_unc)*(n_unc+F*[dq1;dq2])
        ];

    g =[0,0;
        0,0;
        inv(M_unc)];

    dx=f+g*[u1;0];

    q1_0  = x_0(1);
    q2_0  = x_0(2);
    dq1_0 = x_0(3);
    dq2_0 = x_0(4);
    
    m1_var = 1.7851;
    m2_var = 1.9531;

    x = zeros(4,time);
    dstate = zeros(4,time);

    u = zeros(2,time);              %the final input
    tau = zeros(1,time);            %LQR output
    
    D = zeros(7,time);              %experiment output
    h = zeros(1,time);
    h_dot_diff = zeros(1,time);
    h_dot_nom = zeros(1,time);
    q_des = zeros(2,time);
    dstate(:,1) = subs(dx,[q1,q2,dq1,dq2,u1,m1_unc,m2_unc],[q1_0,q2_0,dq1_0,dq2_0,0,m1_var,m2_var]);
    x(:,1) = [q1_0,q2_0,dq1_0,dq2_0]';
    
    if(show_animation == 1)
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
        drawnow;
    end

    for i = 1:time
        q_des(:,i) = Q_DES;

        tau(i) = K_p*(x(1:2,i)-Q_DES) + K_d*x(3:4,i) + tau_eq;
        
        if num_exp == 1 %first experiment executed with k0, CBF-QP controller
            [h(i), h_dot_diff(i), h_dot_nom(i), u(1,i)] = CBFcontroller(x(:,i),tau(i));
        else
            [h(i), h_dot_diff(i), h_dot_nom(i), u(1,i)] =  LCBFcontroller(x(:,i),tau(i),a_est_net,b_est_net);
        end
        
        dstate(:,i) = subs(dx,[q1,q2,dq1,dq2,u1,m1_unc,m2_unc],[x(:,i)',u(1,i),m1_var,m2_var]);
        x(:,i+1) = x(:,i) + dt*dstate(:,i);
        
        D(1:4,i) = x(:,i);
        D(5,i) = u(1,i);
        D(6,i) = h_dot_nom(i);
        D(7,i) = h_dot_diff(i);
        
        if show_animation == 1
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
    end
    cbf_values(exp_num,:) = h;
    %% plot
    % show the state of the robot
    if show_plots
        N = (1 : time) * dt;
        figure();
        subplot(2,2,1);
        plot(N,x(2,1:time),N,q_des(2,1:time))
        % yline(Q2_MAX,'--');
        % yline(-Q2_MAX,'--');
        xlabel('s')
        ylabel('q2')
        legend('q2','q2_des');
        title('q2 angle')
        ylim([-0.5 0.4])


        subplot(2,2,2);
        plot(N,x(1,1:time),N,q_des(1,1:time))
        % yline(pi+Q1_MAX);
        % yline(pi-Q1_MAX);
        ylim([2.5 3.4])
        xlabel('s')
        ylabel('q1')
        legend('q1','q1_des');
        title('q1 angle')

        subplot(2,2,3)
        plot(N,x(3,1:time),N,x(4,1:time));
        yline(1);
        yline(-1);
        xlabel('time');
        ylabel('dq');
        legend('dq1','dq2');
        title('dq velocity');
        ylim([-1.5 1.5])

        subplot(2,2,4)
        plot(x(2,1:time),x(4,1:time))
        yline(1);
        yline(-1);
        xlabel('q2')
        ylabel('dq2')
        title('q2 and dq2')
        ylim([-1.5 1.5])
        grid on

        figure;
        subplot(1,2,1)
        plot(N,u(1,:),N,tau)
        xlabel('time')
        ylabel('input torque')
        legend('u','LQR')
        title('Final control input vs LQR component')

        subplot(1,2,2)
        plot(N,h)
        xlabel('time')
        ylabel('h')
        title('Control Barrier function')
        grid on
    end
end

function [h, h_dot_diff, h_dot_nom, u] = CBFcontroller(x,u_des)
    global a1 a2 a3 a4 a5 f1 f2 dt h_prev
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
    cbf_dot = -c_param*dq2;
    grad = [0,0,0,-c_param*dq2];

    H = 1;
    f_qp = -u_des;

    A = -grad*g;
    b = grad*f+alpha*cbf;

    % disp("A:"); disp(A);
    % disp("b:"); disp(b);
    options = optimset('display','off');
    h = cbf;
    h_dot_diff = (h - h_prev)/dt;
    h_dot_nom = cbf_dot;
    h_prev = h;
    u = quadprog(H,f_qp,A,b,[],[],[],[],[],options);
end

function [h,h_dot_diff, h_dot_nom ,u] = LCBFcontroller(x,u_des,a_est_net,b_est_net)
    global a1 a2 a3 a4 a5 f1 f2 dt h_prev
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
    cbf_dot = -c_param*dq2;
    grad = [0,0,0,-c_param*dq2];

    H = 1;
    f_qp = -u_des;

    a_est = double(predict(a_est_net,x'));
    b_est = double(predict(b_est_net,x'));

    A = -(grad*g + a_est);
    b = grad*f + alpha*cbf + b_est;

    % disp("A:"); disp(A);
    % disp("b:"); disp(b);
    options = optimset('display','off');
    h = cbf;
    h_dot_diff = (h - h_prev)/dt;
    h_dot_nom = cbf_dot;
    h_prev = h;
    u = quadprog(H,f_qp,A,b,[],[],[],[],[],options);
end

function [a_est_net, b_est_net] = erm(D)
    trainingDatastore = getTrainingDataset(D);

    num_state_features = 4;
    num_input_features = 1;
    num_h_dot_features = 1;
    num_hidden_layers = 200;
    input_layer_state = featureInputLayer(num_state_features, 'Name', 'input_state');
    input_layer_u = featureInputLayer(num_input_features, 'Name', 'input_u');
    input_layer_h_dot_nom = featureInputLayer(num_h_dot_features, 'Name', 'input_h_dot_nom');

    layers_a = [
        fullyConnectedLayer(num_hidden_layers, 'Name', 'full_1a')
        tanhLayer('Name','tanh_a')
        fullyConnectedLayer(1,'Name', 'full_2a')];

    layers_b = [
        fullyConnectedLayer(num_hidden_layers, 'Name', 'full_1b')
        tanhLayer('Name','tanh_b')
        fullyConnectedLayer(1,'Name', 'full_2b')];

    mult = multiplicationLayer(2,'Name', 'mult');
    add = additionLayer(3, 'Name', 'add');

    output_layer = regressionLayer('Name','routput');

    lgraph = layerGraph;
    lgraph = addLayers(lgraph, input_layer_state);
    lgraph = addLayers(lgraph, input_layer_u);
    lgraph = addLayers(lgraph, input_layer_h_dot_nom);
    lgraph = addLayers(lgraph,layers_a);
    lgraph = addLayers(lgraph,layers_b);
    lgraph = addLayers(lgraph,mult);
    lgraph = addLayers(lgraph,add);
    lgraph = addLayers(lgraph,output_layer);

    lgraph = connectLayers(lgraph,'input_state','full_1a');
    lgraph = connectLayers(lgraph,'input_state','full_1b');
    lgraph = connectLayers(lgraph,'input_u','mult/in1');
    lgraph = connectLayers(lgraph,'full_2a','mult/in2');
    lgraph = connectLayers(lgraph,'full_2b','add/in1');
    lgraph = connectLayers(lgraph,'mult/out','add/in2');
    lgraph = connectLayers(lgraph,'input_h_dot_nom','add/in3');
    lgraph = connectLayers(lgraph,'add/out','routput');

%     figure
%     plot(lgraph);
%     analyzeNetwork(lgraph);

    % Training settings and train itself
    maxEpochs = 3000;
    miniBatchSize = 1000;
    learningRate = 1e-03;

    options = trainingOptions('adam', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',1, ...
        'InitialLearnRate', learningRate, ...
        'Verbose',true, ...
        'Plots','training-progress');

    net = trainNetwork(trainingDatastore,lgraph,options);

    layers_a = [
        net.Layers(1) %Feature input layer
        net.Layers(4) %First fully connected layer
        net.Layers(5) %Tanh layer
        net.Layers(6) %Second fully connected layer
        net.Layers(12)]; %Output layer

    layers_b = [
        net.Layers(1) %Feature input layer
        net.Layers(7) %First fully connected layer
        net.Layers(8) %Tanh layer
        net.Layers(9) %Second fully connected layer
        net.Layers(12)]; %Output layer

    a_est_net = assembleNetwork(layers_a);
    b_est_net = assembleNetwork(layers_b);
end

function trainingDatastore = getTrainingDataset(D)
    state = D(1:4,:);
    input_u = D(5,:);
    h_dot_nom = D(6,:);
    output = D(7,:);
    
    dataset_dim = size(D);
    
    stateCells = mat2cell(state,4,ones(dataset_dim(2),1));
    inputUCells = mat2cell(input_u,1,ones(dataset_dim(2),1));
    hDotNomCells = mat2cell(h_dot_nom,1,ones(dataset_dim(2),1));
    outputCells = mat2cell(output,1,ones(dataset_dim(2),1));
    stateCells = reshape(stateCells, [dataset_dim(2) 1 1]);
    inputUCells = reshape(inputUCells, [dataset_dim(2) 1 1]);
    hDotNomCells = reshape(hDotNomCells, [dataset_dim(2) 1 1]);
    outputCells = reshape(outputCells, [dataset_dim(2) 1 1]);
    combinedCells = [stateCells inputUCells hDotNomCells outputCells];

    save('trainingData.mat','combinedCells');
    filedatastore = fileDatastore('trainingData.mat','ReadFcn',@load);
    trainingDatastore = transform(filedatastore,@rearrangeData);
end

function out = rearrangeData(ds)
    out = ds.combinedCells;
end