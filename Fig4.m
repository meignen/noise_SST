clc;close all; clear;
load('signal_noise.mat','sn');  
N  = length(sn);
t  = (0:N-1)/N;

%% Computes FSST, FSST2 and FSST3
gamma=0;
sigma = 0.03;
ft =1:N/2;
bt=1:N;
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega2n_1,omega2n_2,omega3,omega3n_1] = sstn_simple(sn,gamma,sigma,ft,bt);
figure;
imagesc(t,ft-1,abs(FSST3));axis xy;xlabel('time','FontSize',20);ylabel('frequency','FontSize',20);title('FSST3','FontSize',20);
set(gca,'ydir','normal')
ax = gca;
ax.FontSize = 20;
figure 
imagesc(t,ft-1,abs(FSST3))
xlim([0.4 0.6]);
ylim([200 300]);
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20)
ax = gca;
ax.FontSize = 20;
hold on; 
[Cs,~] = exridge_mult(abs(FSST2),1,0,10);
plot(t,(Cs(1,:)-1),'LineWidth',4)
plot(t,100+300*t,'--','Linewidth',4);
[Cs,~] = exridge_mult(abs(FSST3),1,0,10);
plot(t,(Cs(1,:)-1),'LineWidth',4)
hold off
