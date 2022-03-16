clc;clear;close all;
load('signal_interference_noise.mat','sn');  
N  = length(sn);
t  = (0:N-1)/N;

% the points bar t_k
tk = (2*(0:floor((4*(260-230)-1)/2))+1)/(4*(260-230));
gamma=0; sigma = 0.03;
ft =1:N/2;bt=1:N;
[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega2n_1,omega2n_2,omega3,omega3n_1] = sstn_simple(sn,gamma,sigma,ft,bt);
xn=0.4;xN=0.5;

[Cs0,] = exridge_mult(abs(STFT), 2, 0, 10);
[Cs,~] = exridge_mult(abs(FSST), 2, 0, 10);
[Cs1,~] = exridge_mult(abs(FSST2), 2, 0, 10);
Z = zeros(1,N);
X = zeros(1,N);
Y = zeros(1,N);
Z1 = zeros(1,N);
X1 = zeros(1,N);
Y1 = zeros(1,N);

for k=1:N
  Z(k)  =  abs(STFT(Cs0(1,k),k));
  X(k)  =  abs(FSST(Cs(1,k),k));
  Y(k)  =  abs(FSST2(Cs1(1,k),k));
  Z1(k) =  abs(STFT(Cs0(2,k),k));
  X1(k) =  abs(FSST(Cs(2,k),k));
  Y1(k) =  abs(FSST2(Cs1(2,k),k));
end

[Loc]  = peakfinder(X);
[Loc1] = peakfinder(Y);
[Loc2] = peakfinder(Z); 

[Loc_1]  = peakfinder(X1);
[Loc1_1] = peakfinder(Y1);
[Loc2_1] = peakfinder(Z1);

ZZ  = zeros(1,length(Loc));
for p = 1:length(Loc)
 ZZ(p) =  abs(FSST(Cs(1,Loc(p)),Loc(p)));
end

ZZ1 = zeros(1,length(Loc1));
for p = 1:length(Loc1)
 ZZ1(p) = abs(FSST2(Cs1(1,Loc1(p)),Loc1(p)));
end

ZZ2 =zeros(1,length(Loc2));
for p = 1:length(Loc2)
 ZZ2(p) = abs(STFT(Cs0(1,Loc2(p)),Loc2(p)));
end

figure
plot((0:N-1)/N,X,(0:N-1)/N,Y,'--',(Loc-1)/N,ZZ,'o',(Loc1-1)/N,ZZ1,'s','Linewidth',2,'MarkerSize',15)
ax = gca;
ax.FontSize = 30;
xlim([xn xN]);
xlabel('time','FontSize',30);
ylabel('amplitude','FontSize',30)
figure
plot((Loc_1-1)/N,Cs(2,Loc_1)-1,'o',(Loc1_1-1)/N,Cs(2,Loc1_1)-1,'s',(Loc2_1-1)/N,Cs0(2,Loc2_1)-1,'*',tk,230*ones(1,length(tk)),'d',(0:N-1)/N,Cs(2,:)-1,(0:N-1)/N,Cs1(2,:)-1,'--',(0:N-1)/N,Cs0(2,:)-1,'-.','Linewidth',2,'MarkerSize',15)
ax = gca;
ax.FontSize = 30;
xlim([xn xN])
ylim([227 235]);
xlabel('time','FontSize',30);
ylabel('frequency','FontSize',30)
figure
plot((Loc-1)/N,Cs(1,Loc)-1,'o',(Loc1-1)/N,Cs(1,Loc1)-1,'s',(Loc2-1)/N,Cs0(1,Loc2)-1,'*',tk,260*ones(1,length(tk)),'d',(0:N-1)/N,Cs(1,:)-1,(0:N-1)/N,Cs1(1,:)-1,'--',(0:N-1)/N,Cs0(1,:)-1,'-.','Linewidth',2,'MarkerSize',15)
ax = gca;
ax.FontSize = 30;
xlim([xn xN])
ylim([255 263]);
xlabel('time','FontSize',30);
ylabel('frequency','FontSize',30)

