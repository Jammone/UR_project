syms x theta dx dtheta u1 real

%% State Variables (Initial Conditions)
x_0 = 4.5;
theta_0 = pi;
dx_0 = 0;
dtheta_0 = 0;

X_EQ = 0;
THETA_EQ = pi;

THETAE_MAX = pi/12;

% Create set of initial conditions

x_init = [x_0,theta_0,dx_0, dtheta_0]'


N = 100
X_0 = zeros(N,4)
%%%% Initialize X_0 with initial conditions for each experiment t=1,..,T
% for i = 1:N
%     X_0(:,i) = x_init
% end

%% Dataset Aggregation for Control Barrier Functions (DaCBarF) -- Episode Learning

%% Inputs

%% CFB (h) {evaluated in x_init}
h = 0.5*(THETAE_MAX^2-(theta-THETA_EQ)^2-c*dtheta^2);

%% Time Derivative Estimate of CBF (h_hat)
grad = [0, THETA_EQ - theta, 0, -c*dtheta];

f =[dx;
    dtheta;
    -inv(M)*(c*[dx;dtheta]+e)
    ];

g =[0,0;
    0,0;
    inv(M)];

h_hat =grad*(f+g*[u1;0]);

%% Nominal State-Feedback Controller
k_0 % from LQR
% This should be the output of LQR (?)

%% Number of Experiments
T = 50

%% Sequence of Truest Coefficients
W = zeros(T)
for i = 2:T-1
    seros(i) = i/T
end


%% Model Classes H_a , H_b

%(OLD VERSION)% Training of Keras NN model in Matlab

% Load Keras model from .json file
% model_a = importKerasLayers('model_a.json')
% model_b = importKerasLayers('model_b.json')

%% (LATEST VERSION) MATLAB DEEP LEARNING TOOL

featuresDim = 4         % Features (input) lenght
hiddenSize = 200;       % Units number of hidden layer
miniBatchSize  = 1      % Mini-batch size

% (inside option -- not necessaty)
%'ValidationData',{XValidation,YValidation},
%'ValidationFrequency',validationFrequency,

options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',30, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20, ...
    'Shuffle','every-epoch', ...
    'Plots','none', ... %'training-progress', ...
    'Verbose',false);

H_a = maeMLP(featuresDim, hiddenSize)
H_b = maeMLP(featuresDim, hiddenSize)

net_a = trainNetwork(XTrain,YTrain,H_a,options);
net_b = trainNetwork(XTrain,YTrain,H_b,options);
output = predict(net, XTrain)



%% Outputs
% Dataset
D
% TIme Derivative Estimate of CBF at last experiment h_hat_T
h_hat_T

% Augmented controller k_t
k_t



function layers = maeMLP(featuresDim, hiddenSize)
    layers = [
        featureInputLayer(numFeatures,'Normalization', 'zscore','name','inputs')
        fullyConnectedLayer(hiddenSize, 'name','fc1')
        batchNormalizationLayer('name','batchNorm')
        reluLayer('name','relu')
        dropoutLayer(0.2,'name','dropout')
        fullyConnectedLayer(1,'name','output')
        maeRegressionLayer('mae')];
end