%% DFT of a sequence and plot of the magnitude and phase response
x = ones(1,4); % Input sequence
N1 = 8; %Number of frequency points
Y1 = dft(x,N1)
k = 0:1:N1-1;
subplot(2,2,1), stem(k,abs(Y1)),   xlabel('k'), ylabel('|Y1(k)|');
subplot(2,2,3), stem(k,angle(Y1)), xlabel('k'), ylabel('arg(Y1(k))');
N2 = 31 %Number of frequency points
Y2 = dft(x,N2)
k = 0:1:N2-1;
subplot(2,2,2), stem(k,abs(Y2)),   xlabel('k'), ylabel('|Y2(k)|');
subplot(2,2,4), stem(k,angle(Y2)), xlabel('k'), ylabel('arg(Y2(k))');

%% Inverse DFT of a sequence
X = [4,1+i,0,1,-i,0,1+i,1-i];
N = length(X);
xn = idft(X,N)

%% Circular convolution of two sequences
n = 0:7;
x = sin(3*pi*n/8); % Input sequence 1
h = [1,1,1,1]; % Input sequence 2
Nx = length(x);
Nh = length(h);
N = 8;
if(N<max(Nx,Nh))
    error('N must be >=max(Nx,Nh')
end
y = circconv(x,h,N)

%% To compare circular and linear convolution of two sequences
x = [1,1,1,2,1,1]; %input sequence
h = [1,1,2,1]; %impulse sequence
Nx = length(x);
Nh = length(h);
N = max(Nx,Nh);
yc = circconv(x,h,N);
y = conv(x,h);
n = 0:1:Nx-1;
subplot(2,2,1), stem(n,x),  xlabel('n'), ylabel('x(n)'),  title('Input Sequence')
n = 0:1:Nh-1;
subplot(2,2,2), stem(n,h),  xlabel('n'), ylabel('h(n)'),  title('Impulse Sequence')
n = 0:1:N-1;
subplot(2,2,3), stem(n,yc), xlabel('n'), ylabel('yc(n)'), title('Circular Convolution')
n = 0:1:Nx+Nh-2;
subplot(2,2,4), stem(n,y),  xlabel('n'), ylabel('y(n)'),  title('Linear Convolution')

%% Overlap and save method
x = [1,2,-1,2,3,-2,-3,-1,1,1,2,-1]; % Input sequence
h = [1,2,1,1]; % Impulse sequence
N = 4; % Length of each block before appending zeros
y = ovrlsav(x,h,N);

%% Overlap and add method
x = [1,2,-1,2,3,-2,-3,-1,1,1,2,-1]; %Input sequence
h = [1,2,1,1]; %Impulse sequence
L = 4; %Length of each block before appending zeros
y = ovrladd(x,h,L);

%% To design a Butterworth lowpass filter for the specifications
alphap = .4;  %Passband attenuation in dB
alphas = 30;  %Stopband attenuation in dB
fp = 400;     %Passband frequency in Hz
fs = 800;     %Stopband frequency in Hz
F = 2000;     %Sampling frequency in Hz
omp = 2*fp/F;
oms = 2*fs/F;
%To find cutoff frequency and order of the filter
[n,wn] = buttord(omp,oms,alphap,alphas);
%system function of the filter
[b,a] = butter(n,wn)
w = 0:.01:pi;
[h,om] = freqz(b,a,w,'whole');
m = abs(h);
an = angle(h);
subplot(2,1,1); plot(om/pi,20*log(m)); grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(om/pi,an);        grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Butterworth bandpass filter for the specifications
alphap = 2;         %Pass band attenuation in dB
alphas = 20;        %Stop band attenuation in dB
wp = [.2*pi,.4*pi]; %Passband frequency in radians
ws = [.1*pi,.5*pi]; %Stopband frequency in radians
%To find cutoff frequency and order of the filter
[n,wn] = buttord(wp/pi,ws/pi,alphap,alphas);
%System function of the filter
[b,a] = butter(n,wn)
w = 0:.01:pi;
[h,ph] = freqz(b,a,w);
m = 20*log10(abs(h));
an = angle(h);
subplot(2,1,1); plot(ph/pi,m); grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(ph/pi,an); grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Butterworth highpass filter for the specifications
alphap = .4;        %Pass band attenuation in dB
alphas = 30;        %Stop band attenuation in dB
fp = 800;           %Passband frequency in radians
fs = 400;           %Stopband frequency in radians
F = 2000;           %Sampling frequency in Hz
omp = 2*fp/F;
oms = 2*fs/F;
%To find cutoff frequency and order of the filter
[n,wn] = buttord(omp,oms,alphap,alphas);
%system function of the filter
[b,a] = butter(n,wn,'high')
w = 0:.01:pi;
[h,om] = freqz(b,a,w);
m = 20*log10(abs(h));
an = angle(h);
subplot(2,1,1); plot(om/pi,m);  grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(om/pi,an); grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Butterworth bandstop filter for the specifications
alphap = 2;         %Pass band attenuation in dB
alphas = 20;        %Stop band attenuation in dB
ws = [.2*pi,.4*pi]; %Stopband frequency in radians
wp = [.1*pi,.5*pi]; %Passband frequency in radians
%To find cutoff frequency and order of the filter
[n,wn] = buttord(wp/pi,ws/pi,alphap,alphas);
%System function of the filter
[b,a] = butter(n,wn,'stop')
w = 0:.01:pi;
[h,ph] = freqz(b,a,w);
m = 20*log10(abs(h));
an = angle(h);
subplot(2,1,1); plot(ph/pi,m);  grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(ph/pi,an); grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Chebyshev 1 lowpass filter for the specifications
alphap = 1;  %Pass band attenuation in dB
alphas = 15; %Stop band attenuation in dB
wp = .2*pi;  %Pass band frequency in radians
ws = .3*pi;  %Stop band frequency in radians
%To find cutoff frequency and order of the filter
[n,wn] = cheb1ord(wp/pi,ws/pi,alphap,alphas);
%System function of the filter
[b,a] = cheby1(n,alphap,wn)
w = 0:.01:pi;
[h,ph] = freqz(b,a,w);
m = 20*log(abs(h));
an = angle(h);
subplot(2,1,1); plot(ph/pi,m);  grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(ph/pi,an); grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Chebyshev 2 lowpass filter for the specifications
alphap = 1;  %Pass band attenuation in dB
alphas = 20; %Stop band attenuation in dB
wp = .2*pi;  %Pass band frequency in radians
ws = .3*pi;  %Stop band frequency in radians
%To find cutoff frequency and order of the filter
[n,wn] = cheb2ord(wp/pi,ws/pi,alphap,alphas);
%System function of the filter
[b,a] = cheby2(n,alphas,wn)
w = 0:.01:pi;
[h,ph] = freqz(b,a,w);
m = abs(h);
an = angle(h);
subplot(2,1,1); plot(ph/pi,20*log(m)); grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(ph/pi,an);        grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Chebyshev 1 bandpass filter for the specifications
alphap = 2;  %Pass band attenuation in dB
alphas = 20; %Stop band attenuation in dB
wp = [.2*pi,.4*pi];  %Pass band frequency in radians
ws = [.1*pi,.5*pi];  %Stop band frequency in radians
%To find cutoff frequency and order of the filter
[n,wn] = cheb1ord(wp/pi,ws/pi,alphap,alphas);
%System function of the filter
[b,a] = cheby1(n,alphap,wn)
w = 0:.01:pi;
[h,ph] = freqz(b,a,w);
m = 20*log10(abs(h));
an = angle(h);
subplot(2,1,1); plot(ph/pi,m);  grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(ph/pi,an); grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To design a Chebyshev 2 bandstop filter for the specifications
alphap = 2;  %Pass band attenuation in dB
alphas = 20; %Stop band attenuation in dB
ws = [.2*pi,.4*pi];  %Stop band frequency in radians
wp = [.1*pi,.5*pi];  %Pass band frequency in radians
%To find cutoff frequency and order of the filter
[n,wn] = cheb2ord(wp/pi,ws/pi,alphap,alphas);
%System function of the filter
[b,a] = cheby2(n,alphas,wn,'stop')
w = 0:.01:pi;
[h,ph] = freqz(b,a,w);
m = 20*log(abs(h));
an = angle(h);
subplot(2,1,1); plot(ph/pi,m);  grid; xlabel('Normalized frequency'); ylabel('Gain in dB');
subplot(2,1,2); plot(ph/pi,an); grid; xlabel('Normalized frequency'); ylabel('Phase in radians');

%% To convert the analog filter into digital filter using impulse invariance
b = [1,2];       % Numerator coefficients of analog filter
a = [1,5,11,15]; % Denominator coefficients of analog filter
f = 5;           % Sampling frequency
[bz,az] = impinvar(b,a,f)

%% To convert the analog filter into digital filter using bilinear transformation
b = [2];     % Numerator coefficients of analog filter
a = [1,3,2]; % Denominator coefficients of analog filter
f = 1;       % Sampling frequency
[bz,az] = bilinear(b,a,f)
