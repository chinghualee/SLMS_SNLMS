clc
clear all
close all

channelType = 1; % 1: quasi // 2: sparse // 3: dispersive

% Load channel impulse response
switch channelType
    case 1
        channelName = 'quasi'
        load('./IRs/ch_quasi.mat')
        channelFilter = ch_quasi;
        N = 1500; % length of signal
    case 2
        channelName = 'sparse'
        load('./IRs/ch_sparse.mat')
        channelFilter = ch_sparse;
        N = 1500; % length of signal
    case 3
        channelName = 'dispersive'
        load('./IRs/ch_dispersive.mat')
        channelFilter = ch_dispersive;
        N = 4000; % length of signal
end
channelFilterLen = length(channelFilter);

% SLMS parameters
adaptiveFilterLen = 256;
mu = 0.0025; % step size parameter
p = [1,1.2,1.5,1.8,2]; % values for p-norm
c = 1e-3; % regularization constant

numTrials = 100; % number of independent Monte Carlo trials

SE = zeros(N,length(p)); % squared error for all cases
for k = 1 : length(p)    
    se = 0; % squared error
    for m = 1 : numTrials

        disp('p ='); disp(p(k)); % display p value
        disp('trial'); disp(m); % display m-th trial
    
        % Generate independent input signals
        u = normrnd(0,1,N,1);
        v = normrnd(0,sqrt(0.01),N,1);

        % Initialize adaptive filter
        adaptiveFilter = zeros(adaptiveFilterLen,1);

        % Buffers for filtering
        channelFilterBuffer = normrnd(0,1,channelFilterLen,1);
        adaptiveFilterBuffer = channelFilterBuffer;

        % u vector for adaptation
        u_vec = zeros(adaptiveFilterLen,1);

        % Sparsity-promoting factors
        s = zeros(adaptiveFilterLen,1);

        % Begin simulation
        for n = 1 : N
            % Obtain the output of channel
            channelFilterBuffer = circshift(channelFilterBuffer,-1,1);
            channelFilterBuffer(end) = u(n);
            y(n) = flip(channelFilter')*channelFilterBuffer;

            % d(n) is the sum of channel output and noise
            d(n) = y(n)+v(n);

            % Obtain output estimate
            adaptiveFilterBuffer = circshift(adaptiveFilterBuffer,-1,1);
            adaptiveFilterBuffer(end) = u(n);
            y_hat(n) = flip(adaptiveFilter')*adaptiveFilterBuffer;

            % e(n) the error signal
            e(n) = d(n)-y_hat(n);

            % Sparsity-promoting factors update
            w = abs(adaptiveFilter)+c;
            s = (2/p(k))*w.^(2-p(k));  
            s = s/mean(s);
            
            % Adaptive filter update
            u_vec = circshift(u_vec,1,1);
            u_vec(1) = u(n);
            adaptiveFilter = adaptiveFilter+mu*s.*u_vec*e(n);
        end

        % Add square error
        se = se+e.^2;
    end

    % Averaging to get MSE
    MSE(:,k) = se/numTrials;
end

MSE = smoothdata(MSE); % smooth the curve for better visualization

figure
hold on
plot(MSE,'lineWidth',2,'lineWidth',2)
set(gca, 'YScale', 'log')
set(gca,'fontSize',16)
set(gcf,'position',[100 100 900 300])
axis([0 N 1e-2 0.2])
title(['SLMS, ' channelName])
xlabel('Iteration ({\it{n}})')
ylabel('MSE')
legend('\it{p}=1.0','\it{p}=1.2','\it{p}=1.5','\it{p}=1.8','\it{p}=2.0')
grid on