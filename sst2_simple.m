function [STFT,SST1,SST2,omega,omega2] = sst2_simple(s,gamma,sigma,ft,bt)
% sstn : computes the STFT of a signal and different versions of synchrosqueezing/reassignment.
%   Uses a Gaussian window, the filter length is the same as the signal
%
% INPUTS:  
%   s: real or complex signal, must be of length 2^N
%   gamma: threshold
%   sigma: window parameter
%   ft: frequency bins
%   bt: time bins
%
% OUTPUTS:   
%   STFT: the short-time Fourier transform
%   SST1: standard synchrosqueezing
%   SST2: vertical second-order synchrosqueezing
%   omega: instantaneous frequency (vertical reassignment operator)
%   omega2: second-order instantaneous frequency

% checking length of signal
n = length(s);
nv = log2(n);
if mod(nv,1)~=0
    warning('The signal is not a power of two, truncation to the next power');
    s = s(1:2^floor(nv));
end
n = length(s);
s = s(:);

% Optional parameters
if nargin<5
   ft = 1:n/2;
   bt = 1:n;
end
nb = length(bt);
neta = length(ft);
sz=zeros(n,1);
% Padding
x =[sz(2:n/2+1); s ; sz(end-n/2:end-1)];

% Window definition
t = -0.5:1/n:0.5-1/n;t=t';
g =  1/sigma*exp(-pi/sigma^2*t.^2);
gp = -2*pi/sigma^2*t .* g; % g'

% Initialization
STFT = zeros(neta,nb);
SST1 = zeros(neta,nb);
SST2 = zeros(neta,nb);
omega = zeros(neta,nb);
tau2 = zeros(neta,nb);
omega2 = zeros(neta,nb);
phi22p = zeros(neta,nb);
vg = zeros(neta,7);
vgp = zeros(neta,5);

%% Computes STFT and reassignment operators
for b=1:nb
	
    % STFT, window t^n*g
    for i = 0:2
        tmp = (fft(x(bt(b):bt(b)+n-1).*(t.^i).*g))/n;
        vg(:,i+1) = tmp(ft);
    end 
       
      %% STFT, window gx^np
    for i = 0:1
        tmp = fft(x(bt(b):bt(b)+n-1).*(t.^i).*gp)/n;
        vgp(:,i+1) = tmp(ft);
    end
    
    %% second-order operator tau
    tau2(:,b) = vg(:,2)./vg(:,1);
    
    %% operator omega
    omega(:,b) = (ft-1)'-real(vgp(:,1)/2/1i/pi./vg(:,1));

    %% operator hat p: estimations of frequency modulation  
    W2 = 1/2/1i/pi*(vg(:,1).^2+vg(:,1).*vgp(:,2)-vg(:,2).*vgp(:,1));
    Y2 = vg(:,1).*vg(:,3) - vg(:,2).*vg(:,2);
    phi22p(:,b) = W2./Y2;
    omega2(:,b) = omega(:,b) + real(phi22p(:,b).*tau2(:,b));
     
    % Storing STFT
    STFT(:,b) = vg(:,1).* exp(1i*pi*(ft-1)'); % compensates the tranlation 1/2 of s
        
end

%% Reassignment step
for b=1:nb
    for eta=1:neta
        if abs(STFT(eta,b))>gamma
           
%%%%%%%%%%%%%%%%%%%%%%%%%%SST1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            k = 1+round(omega(eta,b));          
            if k>=1 && k<=neta               
             SST1(k,b)  = SST1(k,b) + STFT(eta,b);
            end
            
%%%%%%%%%%%%%%%%%%%%SST2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            k = 1+round(omega2(eta,b));
             if k>=1 && k<=neta
               SST2(k,b)  = SST2(k,b) + STFT(eta,b);
             end
             
        end
    end
end
