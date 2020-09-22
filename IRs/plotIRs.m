% Plot impulse responses (IRs)

load('ch_quasi.mat')
figure
plot(ch_quasi)
title('Quasi-sparse')
xlabel('Time (samples)')
ylabel('Amplitude')
axis([0 256 -0.18 0.18])
grid on
set(gca,'fontSize',16)
set(gcf,'position',[500 500 900 300])

load('ch_sparse.mat')
figure
plot(ch_sparse)
title('Sparse')
xlabel('Time (samples)')
ylabel('Amplitude')
axis([0 256 -0.3 0.3])
grid on
set(gca,'fontSize',16)
set(gcf,'position',[500 500 900 300])

load('ch_dispersive.mat')
figure
plot(ch_dispersive)
title('Dispersive')
xlabel('Time (samples)')
ylabel('Amplitude')
axis([0 256 -0.08 0.08])
grid on
set(gca,'fontSize',16)
set(gcf,'position',[500 500 900 300])

