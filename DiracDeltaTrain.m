function [y]  = DiracDeltaTrain(duration,period,numExp)

% duration -> duration of the impulse train in seconds
% period -> period of the impulse train in seconds
% numExp -> Fourier series representation of periodic DiracDelta train 
% is used for impulse train generation. This argument is the number of
% exponential terms that are to be summed to produce the train.
% Higher this value is, better the approximation. 

Fs = 10; % sampling frequency in Hz
T = -duration/2 : 1/Fs : duration/2;

sum = zeros(1,length(T));

for i = 1:length(T)
    for k = -numExp : numExp
        sum(i) = sum(i) + (1/period)*exp(-j*2*(pi/period)*k*T(i));
    end
end

plot(T,abs(sum));
title("Periodic Dirac-Delta Impulse train");
xlabel("Time in seconds");ylabel("Magnitude");
