clc;clear;close all;
%% signal
N  = 1024;t  = (0:N-1)/N;fs = 0:N/2;ft =1:N/2;bt=1:N;
phi1 = 260*t;phi2 = 230*t;
s = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
%% variables to compute TFRs
gamma=0;sigma = 0.03;
%% x y limits to zoom the TFRs
xn=0.4;xN=0.6;yn1=180;yN1=300;
%% tilde tk 
k=0:2:length(t);
tz=(k+1)/(2*(260-230));
p=interp1(t, t, tz(20), 'nearest');t_k1=find(t==p);t(t_k1);
%% tk
tz=(k)/(2*(260-230));
p=interp1(t, t, tz(20), 'nearest');t_k2=find(t==p);t(t_k2);
%% Computes TFRs
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega2n_1,omega2n_2,omega3,omega3n_1] = sstn_simple(s,gamma,sigma,ft,bt);
[Cs, Es] = exridge_mult(abs(STFT), 2, 0, 10);

%computation of reassignment vectors at tk
figure
tt  = 20;
tt1 = 10;
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,omega(Cs(2,t_k2)-tt:Cs(2,t_k2)+tt1,t_k2),'-*','Linewidth',2,'MarkerSize',10);
hold on;
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,omega2(Cs(2,t_k2)-tt:Cs(2,t_k2)+tt1,t_k2),'Linewidth',2,'MarkerSize',10);
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,(Cs(2,t_k2)-1)*ones(1,tt+tt1+1),'--','Linewidth',2);
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,230*ones(1,tt+tt1+1),'Linewidth',2);
xlabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
ylim([225 235]);
yticks(225:235);
xlim([Cs(2,t_k2)-1-tt Cs(2,t_k2)-1+5]);
legend({'$\widehat{\omega}_{f}$','$\widehat{\omega}_{f}^{[2]}$','spectrogram ridge','$\xi_1$',...
},'Interpreter','latex');
hold off;

%computation of reassignment vector at tilde tk
figure
tt  = 20;
tt1 = 20;
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,omega(Cs(2,t_k1)-tt:Cs(2,t_k1)+tt1,t_k1),'-*','Linewidth',2,'MarkerSize',10);
hold on;
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,omega2(Cs(2,t_k1)-tt:Cs(2,t_k1)+tt1,t_k1),'Linewidth',2,'MarkerSize',10);
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,(Cs(2,t_k1)-1)*ones(1,tt+tt1+1),'--','Linewidth',2);
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,230*ones(1,tt+tt1+1),'Linewidth',2);
xlabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
ylim([225 235]);
yticks(225:235);
xlim([210 240]);
legend({'$\widehat{\omega}_{f}$','$\widehat{\omega}_{f}^{[2]}$','spectrogram ridge','$\xi_1$',...
},'Interpreter','latex');
hold off;

%computation of second order reassignment vectors at tk
figure
tt  = 10;
tt1 = 10;
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,omega2(Cs(2,t_k2)-tt:Cs(2,t_k2)+tt1,t_k2),'-*','Linewidth',2,'MarkerSize',10);
hold on;
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,omega2n_1(Cs(2,t_k2)-tt:Cs(2,t_k2)+tt1,t_k2),'-o','Linewidth',2,'MarkerSize',10);
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,omega2n_2(Cs(2,t_k2)-tt:Cs(2,t_k2)+tt1,t_k2),'-+','Linewidth',2,'MarkerSize',10);
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,(Cs(2,t_k2)-1)*ones(1,tt+tt1+1),'--','Linewidth',2);
plot(Cs(2,t_k2)-1-tt:Cs(2,t_k2)-1+tt1,230*ones(1,tt+tt1+1),'Linewidth',2);
xlabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
ylim([225 235]);
yticks(225:238);
xlim([225 238]);
legend({'$\widehat{\omega}_{f}^{[2]}$','$\widehat{\omega}_{f,1}^{[2]}$','$\widehat{\omega}_{f,2}^{[2]}$','spectrogram ridge','$\xi_1$',...
},'Interpreter','latex');
hold off;

%computation of reassignment vector at tilde tk
figure
tt  = 10;
tt1 = 5;
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,omega2(Cs(2,t_k1)-tt:Cs(2,t_k1)+tt1,t_k1),'-*','Linewidth',2,'MarkerSize',10);
hold on;
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,omega2n_1(Cs(2,t_k1)-tt:Cs(2,t_k1)+tt1,t_k1),'-o','Linewidth',2,'MarkerSize',10);
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,omega2n_2(Cs(2,t_k1)-tt:Cs(2,t_k1)+tt1,t_k1),'-+','Linewidth',2,'MarkerSize',10);
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,(Cs(2,t_k1)-1)*ones(1,tt+tt1+1),'--','Linewidth',2);
plot(Cs(2,t_k1)-1-tt:Cs(2,t_k1)-1+tt1,230*ones(1,tt+tt1+1),'Linewidth',2);
xlabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
ylim([225 235]);
yticks(225:235);
xlim([Cs(2,t_k1)-1-tt Cs(2,t_k1)-1+tt1])
legend({'$\widehat{\omega}_{f}^{[2]}$','$\widehat{\omega}_{f,1}^{[2]}$','$\widehat{\omega}_{f,2}^{[2]}$','spectrogram ridge','$\xi_1$',...
},'Interpreter','latex');
hold off;
