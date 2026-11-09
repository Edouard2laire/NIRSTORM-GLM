function A = acf_rows(X, maxlag)          % ACF of each row, via FFT
    X = X - mean(X, 2);
    n = size(X, 2);
    nfft = 2^nextpow2(2*n - 1);
    F = fft(X, nfft, 2);
    c = real(ifft(F .* conj(F), nfft, 2));
    A = c(:, 1:maxlag+1) ./ c(:, 1);
end
