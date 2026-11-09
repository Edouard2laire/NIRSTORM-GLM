function [acf, lags] = acf_base(x, maxlag)
    x = x(:) - mean(x);
    n = numel(x);
    nfft = 2^nextpow2(2*n - 1);
    F = fft(x, nfft);
    c = real(ifft(F .* conj(F)));       % autocovariance via FFT
    acf = c(1:maxlag+1) / c(1);         % normalize so lag 0 = 1
    lags = (0:maxlag)';
end