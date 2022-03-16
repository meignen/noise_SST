clc;close all; clear;

%% signal linear chirp 

N  = 1024;
t  = (0:N-1)/N;
fs = 0:N/2;

phi1 = 260*t;
phi2 = 230*t;

s0 = exp(2*pi*1i*(phi1))+exp(2*pi*1i*(phi2));
phi10_prim = 230;
phi20_prim = 260;

nbreal = 40; 
SNR = -5:5;
SNR_IF1   = zeros(1,nbreal); %computation of IF with omega on the ridge of FSST
SNR_IF2   = zeros(1,nbreal); %computation of IF with omega2 on the ridge of FSST2
SNR_IF1_1 = zeros(1,nbreal); %computation of IF with spline approximation using maxima on FSST ridge
SNR_IF1_2 = zeros(1,nbreal); %computation of IF with spline approximation using maxima on FSST ridge
SNR_IF1_3 = zeros(1,nbreal); %computation of IF with spline approximation using maxima on FSST ridge
SNR_IF2_1 = zeros(1,nbreal); %computation of IF with spline approximation using maxima on FSST ridge
SNR_IF2_2 = zeros(1,nbreal); %computation of IF with spline approximation using maxima on FSST ridge
SNR_IF2_3 = zeros(1,nbreal); %computation of IF with spline approximation using maxima on FSST ridge

figure
gamma=0;
sigma = 0.03;
ft =1:N/2;
bt=1:N;
%[STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega3,tau2,tau3,phi22p,phi33p,phi44p] = sstn(s0,gamma,sigma,ft,bt);
[STFT,FSST,FSST2,omega,omega2] = sst2_simple(s0,gamma,sigma,ft,bt);
imagesc(t,ft,abs(STFT))
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;

%we compute the IF approximations for s0
SNRIF1 = zeros(1,length(SNR));
SNRIF2 = zeros(1,length(SNR));
SNRIF11 = zeros(1,length(SNR));
SNRIF12 = zeros(1,length(SNR));
SNRIF13 = zeros(1,length(SNR));
SNRIF21 = zeros(1,length(SNR));
SNRIF22 = zeros(1,length(SNR));
SNRIF23 = zeros(1,length(SNR));


for k=1:length(SNR)
 
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

  %we compute the estimate using omega
  val1    = csaps((Loc-1)/N,omega_r(Loc),1,t,ones(1,length(Loc))); 
  val095  = csaps((Loc-1)/N,omega_r(Loc),0.95,t,ones(1,length(Loc)));
  val09   = csaps((Loc-1)/N,omega_r(Loc),0.9,t,ones(1,length(Loc)));

  %we compute the estimate using omega2
  val1_2    = csaps((Loc1-1)/N,omega2_r(Loc1),1,t,ones(1,length(Loc1))); 
  val095_2  = csaps((Loc1-1)/N,omega2_r(Loc1),0.95,t,ones(1,length(Loc1)));
  val09_2   = csaps((Loc1-1)/N,omega2_r(Loc1),0.9,t,ones(1,length(Loc1)));

  SNR_IF1(p)   = snr(phi10_prim,phi10_prim-omega_r);
  SNR_IF2(p)   = snr(phi10_prim,phi10_prim-omega2_r);
  SNR_IF1_1(p) = snr(phi10_prim,phi10_prim-val1);
  SNR_IF1_2(p) = snr(phi10_prim,phi10_prim-val095);
  SNR_IF1_3(p) = snr(phi10_prim,phi10_prim-val09);
  SNR_IF2_1(p) = snr(phi10_prim,phi10_prim-val1_2);
  SNR_IF2_2(p) = snr(phi10_prim,phi10_prim-val095_2);
  SNR_IF2_3(p) = snr(phi10_prim,phi10_prim-val09_2);
 
 end
 SNRIF1(k)  = mean(SNR_IF1);
 SNRIF2(k)  = mean(SNR_IF2);
 SNRIF11(k) = mean(SNR_IF1_1);
 SNRIF12(k) = mean(SNR_IF1_2);
 SNRIF13(k) = mean(SNR_IF1_3);
 SNRIF21(k) = mean(SNR_IF2_1);
 SNRIF22(k) = mean(SNR_IF2_2);
 SNRIF23(k) = mean(SNR_IF2_3);
end

figure
hold on
plot(SNR,SNRIF1,SNR,SNRIF11,'--',SNR,SNRIF12,'-.',SNR,SNRIF13,':','Linewidth',2);
plot(SNR,SNRIF2,'-d',SNR,SNRIF21,'-s',SNR,SNRIF22,'-o',SNR,SNRIF23,'->','Linewidth',2,'MarkerSize',10);
xlabel('SNR','FontSize',20);
ylabel('SNR out','FontSize',20);
legend({'FSST ridge ,$\widehat{\omega}_{\tilde f}$','FSST,$\hat s$,p=1','FSST,$\hat s$,p=0.95','FSST,$\hat s$,p=0.9',...
    'FSST2 ridge,$\widehat{\omega}_{\tilde f}^{[2]}$','FSST2,$\hat s$,p=1','FSST2,$\hat s$,p=0.95','FSST2,$\hat s$,p=0.9'},...
    'Interpreter','latex');
ax = gca;
ax.FontSize = 20;
hold off
