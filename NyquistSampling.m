% Illustrates Nyquist Sampling Theorem
% =====================================
% 1. Generate a cosine signal of a given frequency.
% 2. Generate the spectrum of this signal.
% 3. Sample the signal at a Nyquist frequency.
% 4. Genrate the DTFT of the sampled signal.
% 5. Under sample the given signal.
% 6. Generate the DTFT of the sampled signal to illustrate aliasing.

% Cosine signal of 8000Hz
Dt = 1/100000;
t = -10/1000:Dt:10/1000;
x = cos(16000*pi*t);

subplot(3,2,1);
plot(t*1000,x);
xlabel("Time in milliseconds");
ylabel("Magintude");
title("Cosine signal of frequency 16000Hz");

hold on;

% Spectrum of the signal -> 'Continuous time' Fourier Transform of the signal
Wmax = 20000 * pi; K = 30000; k = -K:1:K; w = (Wmax/K)*k;
X = x * exp(-j*t'*w) * Dt;

subplot(3,2,2);
plot(w/(2*pi),abs(X));
xlabel("Frequency in Hz");
ylabel("Magintude");
title("Fourier Transform of the Cosine signal");

% Sample the signal at frequency of 20000Hz > Nyquist frequency of 16000Hz
Ts = 1/20000;
n = -100:1:100; 
xs = cos(16000*pi*(n*Ts));

subplot(3,2,3);
stem(n*Ts,xs);
xlabel("Sample Index");
ylabel("Magintude");
title("Sampled Cosine signal (sampled at 20000Hz)");

% DTFT of the sampled signal
Wmax = 2*pi; K = 500; k = -K:1:K; w = (Wmax/K)*k;
X = xs * exp(-j*n'*w);

subplot(3,2,4);
plot(w/pi,abs(X));
xlabel("Digital angular frequency");
ylabel("Magintude");
title("DTFT of Sampled Cosine signal(sampled at 20000Hz)");

% Under-sample the signal at frequency of 6000Hz < Nyquist frequency of 16000Hz
Ts = 1/6000;
n = -60:1:60; 
xus = cos(4000*pi*(n*Ts));

subplot(3,2,5);
stem(n*Ts,xus);
xlabel("Sample Index");
ylabel("Magintude");
title("Under-sampled Cosine signal (sampled at 6000Hz)");

% DTFT of the under-sampled signal
Wmax = 2*pi; K = 500; k = -K:1:K; w = (Wmax/K)*k;
Xu = xus * exp(-j*n'*w);

subplot(3,2,6);
plot(w/pi,abs(Xu));
xlabel("Digital angular frequency");
ylabel("Magintude");
title("DTFT of under-sampled Cosine signal(sampled at 6000Hz)");



