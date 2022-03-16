clc;close all; clear;

load -ascii batsig.txt
s = batsig;
s = s(145:end)';
s = hilbert(s);
s = s(:);
N = length(s);
t  = (0:N-1)/N;
if (N == 2^(floor(log2(N))))
 Nfft = N;
else
 Nfft = 2^(floor(log2(N))+1);
end
sigma_opt = 0.13;  
clwin = 10;
gamma = 10^(-6);
[STFT,FSST,FSST2,omega,omega2] = sst2(s,sigma_opt,Nfft,gamma);
figure
imagesc(t,0:Nfft/2-1,abs(STFT(1:Nfft/2,:)))
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
hold on; 
[Cs, ~] = exridge_mult(abs(FSST),4,0,10);

for k=1:4
  plot(t,Cs(k,:),'Linewidth',2,'MarkerSize',10);
end
hold off;

figure
imagesc(t,0:Nfft/2-1,abs(STFT(1:Nfft/2,:)))
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
hold on; 
[Cs1, Es] = exridge_mult(abs(FSST2),4,0,10);
for k=1:4
  plot(t,Cs1(k,:),'Linewidth',2,'MarkerSize',10);
end
hold off;

X1 = zeros(1,N);
X2 = zeros(1,N);
X3 = zeros(1,N);
X4 = zeros(1,N);
omega2_r1 = zeros(1,N);
omega2_r2 = zeros(1,N);
omega2_r3 = zeros(1,N);
omega2_r4 = zeros(1,N);

for q=1:N
 X1(q) = abs(FSST2(Cs1(1,q),q));
 X2(q) = abs(FSST2(Cs1(2,q),q));
 X3(q) = abs(FSST2(Cs1(3,q),q));
 X4(q) = abs(FSST2(Cs1(4,q),q));
 omega2_r1(q) = omega2(Cs1(1,q),q);
 omega2_r2(q) = omega2(Cs1(2,q),q);
 omega2_r3(q) = omega2(Cs1(3,q),q);
 omega2_r4(q) = omega2(Cs1(4,q),q);
end
 
[Loc]  = peakfinder(X1);
[Loc1] = peakfinder(X2);
[Loc2] = peakfinder(X3);
[Loc3] = peakfinder(X4);

%first mode with omega2 at local mxima
val1_1    = csaps((Loc-1)/N,omega2_r1(Loc),1,t,ones(1,length(Loc))); 
val095_1  = csaps((Loc-1)/N,omega2_r1(Loc),0.95,t,ones(1,length(Loc)));
val09_1   = csaps((Loc-1)/N,omega2_r1(Loc),0.9,t,ones(1,length(Loc)));

%second mode with omega2 at local maxima
val1_2    = csaps((Loc1-1)/N,omega2_r2(Loc1),1,t,ones(1,length(Loc1))); 
val095_2  = csaps((Loc1-1)/N,omega2_r2(Loc1),0.95,t,ones(1,length(Loc1)));
val09_2   = csaps((Loc1-1)/N,omega2_r2(Loc1),0.9,t,ones(1,length(Loc1)));

%third mode with omega2 at local maxima
val1_3    = csaps((Loc2-1)/N,omega2_r3(Loc2),1,t,ones(1,length(Loc2))); 
val095_3  = csaps((Loc2-1)/N,omega2_r3(Loc2),0.95,t,ones(1,length(Loc2)));
val09_3   = csaps((Loc2-1)/N,omega2_r3(Loc2),0.9,t,ones(1,length(Loc2)));

%third mode with omega2 at local maxima
val1_4    = csaps((Loc3-1)/N,omega2_r4(Loc3),1,t,ones(1,length(Loc3))); 
val095_4  = csaps((Loc3-1)/N,omega2_r4(Loc3),0.95,t,ones(1,length(Loc3)));
val09_4   = csaps((Loc3-1)/N,omega2_r4(Loc3),0.9,t,ones(1,length(Loc3)));

figure
imagesc(t,0:Nfft/2-1,abs(STFT(1:Nfft/2,:)))
set(gca,'ydir','normal');
xlabel('time','FontSize',20);
ylabel('frequency','FontSize',20);
ax = gca;
ax.FontSize = 20;
hold on; 
plot(t,val095_1,'Linewidth',2);
plot(t,val095_2,'Linewidth',2);
plot(t,val095_3,'Linewidth',2);
plot(t,val095_4,'Linewidth',2);
plot((Loc-1)/N,omega2_r1(Loc),'d','Markersize',10);
plot((Loc1-1)/N,omega2_r2(Loc1),'d','Markersize',10);
plot((Loc2-1)/N,omega2_r3(Loc2),'d','Markersize',10);
plot((Loc3-1)/N,omega2_r4(Loc3),'d','Markersize',10);


%t,val095_1,'--',t,val09_1,'-.')
hold off
