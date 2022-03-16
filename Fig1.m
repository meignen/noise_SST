clc;clear;close all;
%% signal
N  = 1024;t  = (0:N-1)/N;fs = 0:N/2;ft =1:N/2;bt=1:N;
phi1 = 260*t;phi2 = 230*t;
s = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
gamma=0; sigma = 0.03;
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega2n_1,omega2n_2,omega3,omega3n_1] = sstn_simple(s,gamma,sigma,ft,bt);
[Cs0,] = exridge_mult(abs(STFT), 2, 0, 10);
[Cs,~] = exridge_mult(abs(FSST), 2, 0, 10);
[Cs1,~] = exridge_mult(abs(FSST2), 2, 0, 10);

xn=0.4;xN=0.5;yn1=180;yN1=300;
figure;
imagesc(t,fs,abs(STFT));
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
title('STFT','FontSize',20);
ax = gca;
ax.FontSize = 20;
xlim([xn xN]);
ylim([yn1 yN1]);
hold on;
plot(t,Cs0(1,:),'Linewidth',2);
plot(t,Cs0(2,:),'Linewidth',2);
hold off

figure 
imagesc(t,fs,abs(FSST));
xlim([xn xN]);
ylim([yn1 yN1]);
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
title('FSST','FontSize',20);
ax = gca;
ax.FontSize = 20;
xlim([xn xN]);
ylim([yn1 yN1]);
hold on;
plot(t,Cs(1,:),'Linewidth',2);
plot(t,Cs(2,:),'Linewidth',2);
hold off

figure;
imagesc(t,fs,abs(FSST2));
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
title('FSST2','FontSize',20);
ax = gca;
ax.FontSize = 20;
xlim([xn xN]);
ylim([yn1 yN1]);
hold on;
plot(t,Cs1(1,:),'Linewidth',2);
plot(t,Cs1(2,:),'Linewidth',2);
hold off
