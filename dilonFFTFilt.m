function y = dilonFFTFilt(b, x, nfft)
% fftfilt uses slow overlap/add implementation of fft filtering, we use 
% just use fft/ifft to perform convolution; uses more memory
% b: filter coefficients
% x: data to be filtered by b
% nfft: number of points to be used in fft
%
% NOTE: do not use this function, but copy its functionality, if performing
% the filtering pushes the limits of the machines RAM as copying the data
% for the return will likely cause swapping

% t0 = clock; %execution time, for debug

% nfft = size(x,1);

X = fft(x,nfft);
B = fft(b,nfft);

Y = bsxfun(@times,B,X);

y = ifft(Y);

% t1 = etime(clock,t0);
% disp(['(status) dilonFFTFilt']);
% disp(sprintf('(status) running time: %0.2f sec', t1));