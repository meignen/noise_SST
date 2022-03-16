clc;clear;close all;

%% signal
N  = 1024;t  = (0:N-1)/N;
ft =1:N/2;bt=1:N;
phi = 100*t+150*t.^2;
phi_prim = 100 + 300*t;
s  = exp(2*pi*1i*(phi));

SNR = 5:20;
nb_real = 20;
err1 = zeros(1,length(SNR));
err2 = zeros(1,length(SNR));
err3 = zeros(1,length(SNR));
err4 = zeros(1,length(SNR));
err11 = zeros(1,length(SNR));
err22 = zeros(1,length(SNR));
err33 = zeros(1,length(SNR));
err44 = zeros(1,length(SNR));

gamma=0; sigma = 0.03;
err1_SNR = zeros(1,nb_real);
err2_SNR = zeros(1,nb_real);
err3_SNR = zeros(1,nb_real);
err4_SNR = zeros(1,nb_real);
err11_SNR = zeros(1,nb_real);
err21_SNR = zeros(1,nb_real);
err31_SNR = zeros(1,nb_real);
err41_SNR = zeros(1,nb_real);

for k = 1:length(SNR) 
 
 

 for p = 1:nb_real  
  n = randn(N,1)+1i*randn(N,1);
  [sn] = sigmerge(s(:),n,SNR(k));
  % the points bar t_k
  [STFT,FSST,FSST2,FSST3,FSST4,omega,omega2,omega2n_1,omega2n_2,omega3,omega3n_1]=sstn_simple(sn,gamma,sigma,ft,bt);
  [Cs,~] = exridge_mult(abs(FSST), 1, 0, 10);
  [Cs1,~] = exridge_mult(abs(FSST2), 1, 0, 10);
  Y = zeros(1,N);
  X = zeros(1,N);

  for kk = 1:N 
   X(kk) = abs(FSST(Cs(1,kk),kk));
   Y(kk) = abs(FSST2(Cs1(1,kk),kk));
  end

  [Loc]  = peakfinder(X);
  [Loc1] = peakfinder(Y);


   ZZ  = zeros(1,N);
   for pp = 1:N
    ZZ(pp) = omega(Cs(1,pp),pp); 
   end

   
   ZZ1 =zeros(1,N);
   for pp = 1:N
    ZZ1(pp) = omega2(Cs1(1,pp),pp);
   end

  err1_SNR(p) = mean(abs(phi_prim(Loc)-ZZ(Loc)));
  err2_SNR(p) = mean(abs(phi_prim-ZZ));
  err3_SNR(p) = mean(abs(phi_prim(Loc1)-ZZ1(Loc1)));
  err4_SNR(p) = mean(abs(phi_prim-ZZ1));
  err11_SNR(p) = mean(abs(phi_prim(Loc)-(Cs(Loc)-1)));
  err21_SNR(p) = mean(abs(phi_prim-(Cs(1,:)-1)));
  err31_SNR(p) = mean(abs(phi_prim(Loc1)-(Cs1(Loc1)-1)));
  err41_SNR(p) = mean(abs(phi_prim-(Cs1(1,:)-1)));
 
 
 end   
 err1(k) = mean(err1_SNR);
 err2(k) = mean(err2_SNR);
 err3(k) = mean(err3_SNR);
 err4(k) = mean(err4_SNR);
 err11(k) = mean(err11_SNR);
 err22(k) = mean(err21_SNR);
 err33(k) = mean(err31_SNR);
 err44(k) = mean(err41_SNR);
end
figure
plot(5:20,err1,5:20,err2,'--',5:20,err3,'-.',5:20,err4,':',5:20,err11,'-d',5:20,err22,'->',...
     5:20,err33,'-<', 5:20,err44,'-o','Linewidth',2,'MarkerSize',15);
xlabel('time','FontSize',30);
ylabel('E_{freq}','FontSize',30)

