clc;close all; clear;

%% signal linear chirp 

N  = 1024;
t  = (0:N-1)/N;
ft = 1:N/2;
bt = 1:N;
gamma=0;
sigma = 0.03;

phi1 = 260*t;
phi2 = 230*t;

s0 = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));

nbreal = 40; 
SNR = -5:5;

max_1   = zeros(1,nbreal); %computation of IF with omega on the ridge of FSST
max_2   = zeros(1,nbreal); %computation of IF with omega2 on the ridge of FSST2

%%%%%%%%%%%%%%%%%%%%we compute the IF approximations for s0

maxs0_1 = zeros(1,length(SNR));
maxs0_2 = zeros(1,length(SNR));

for k=1:length(SNR)
 k
 for p = 1:nbreal
  n    = randn(N,1)+1i*randn(N,1);
  [sn]  = sigmerge(s0(:),n,SNR(k));
  
  [STFT,FSST,FSST2,omega,omega2] = sst2_simple(sn,gamma,sigma,ft,bt);

  [Cs,~] = exridge_mult(abs(FSST),2,0,10);
  [Cs1,~] = exridge_mult(abs(FSST2),2,0,10);
  X = zeros(1,N);
  Y = zeros(1,N);
  omega_r  = zeros(1,N);
  omega2_r = zeros(1,N);

  for q=1:N 
   X(q) = abs(FSST(Cs(2,q),q));
   Y(q) = abs(FSST2(Cs1(2,q),q));
   omega_r(q) = omega(Cs(2,q),q);
   omega2_r(q) = omega2(Cs1(2,q),q);
  end

  [Loc]  = peakfinder(X);
  [Loc1] = peakfinder(Y);
  max_1(p) = length(Loc)/N;
  max_2(p) = length(Loc1)/N;  
 end
 
 maxs0_1(k)  = mean(max_1);
 maxs0_2(k)  = mean(max_2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi = 100*t+150*t.^2;
s  = exp(2*pi*1i*(phi));

maxs_1 = zeros(1,length(SNR));
maxs_2 = zeros(1,length(SNR));

%we compute the IF approximations for s0
for k=1:length(SNR)
 k
 for p = 1:nbreal
    n    = randn(N,1)+1i*randn(N,1);
    [sn]  = sigmerge(s(:),n,SNR(k));
  
   %% Computes FSST, FSST2 and FSST3
  [STFT,FSST,FSST2,omega,omega2] = sst2_simple(sn,gamma,sigma,ft,bt);

  [Cs,~]  = exridge_mult(abs(FSST),1,0,10);
  [Cs1,~] = exridge_mult(abs(FSST2),1,0,10);
  
  X = zeros(1,N);
  Y = zeros(1,N);
  omega_r  = zeros(1,N);
  omega2_r = zeros(1,N);

  for q=1:N 
   X(q) = abs(FSST(Cs(1,q),q));
   Y(q) = abs(FSST2(Cs1(1,q),q));
   omega_r(q)  = omega(Cs(1,q),q);
   omega2_r(q) = omega2(Cs1(1,q),q);
  end

  [Loc]  = peakfinder(X);
  [Loc1] = peakfinder(Y);
  
  max_1(p) = length(Loc)/N;
  max_2(p) = length(Loc1)/N;
 end

 maxs_1(k) = mean(max_1);
 maxs_2(k) = mean(max_2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = exp(2*pi*1i*(250*t+20*cos(3*pi*t)));

maxs1_1   = zeros(1,length(SNR)); %computation of IF with omega on the ridge of FSST
maxs1_2   = zeros(1,length(SNR)); %computation of IF with omega2 on the ridge of FSST2

%we compute the IF approximations for s0
for k=1:length(SNR)
 k
 for p = 1:nbreal
  n    = randn(N,1)+1i*randn(N,1);
  [sn]  = sigmerge(s1(:),n,SNR(k));
  
  [STFT,FSST,FSST2,omega,omega2] = sst2_simple(sn,gamma,sigma,ft,bt);

  [Cs,~] = exridge_mult(abs(FSST),1,0,10);
  [Cs1,~] = exridge_mult(abs(FSST2),1,0,10);
  
  X = zeros(1,N);
  Y = zeros(1,N);
  omega_r = zeros(1,N);
  omega2_r = zeros(1,N);

  for q=1:N 
   X(q) = abs(FSST(Cs(1,q),q));
   Y(q) = abs(FSST2(Cs1(1,q),q));
   omega_r(q)  = omega(Cs(1,q),q);
   omega2_r(q) = omega2(Cs1(1,q),q);
  end

  [Loc]  = peakfinder(X);
  [Loc1] = peakfinder(Y);
  max_1(p) = length(Loc)/N;
  max_2(p) = length(Loc1)/N;
 
 end
 maxs1_1(k)  = mean(max_1);
 maxs1_2(k)  = mean(max_2);
end

figure
plot(SNR,maxs0_1,SNR,maxs_1,'--',SNR,maxs1_1,'-.',...
SNR,maxs0_2,'-d',SNR,maxs_2,'-s',SNR,maxs1_2,'-o','Linewidth',2,'MarkerSize',10);
xlabel('SNR','FontSize',20);
ylabel('proportion of modulus maxima','FontSize',20);
ax = gca;
ax.FontSize = 20;

