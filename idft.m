function xn = idft(X,N)
%To compute the inverse DFT of the sequence X(k)
L = length(X); %Length of the sequence
%Computation of twiddle factors
for k = 0:1:N-1
    for n = 0:1:N-1
        p = exp(i*2*pi*n*k/N);
        x2(k+1,n+1) = p;
    end
end
xn = (X*x2.')/N;
