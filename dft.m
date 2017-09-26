function X = dft(xn,N)
% To compute the DFT of the sequence x(n)
L = length(xn); %Length of the sequence
%Checking for the length of the DFT
if(N<L)
    error('N must be >=L')
end
x1 = [xn zeros(1,N-L)]; %Appending zeros
%Computation of twiddle factors
for k = 0:1:N-1
    for n = 0:1:N-1
        p = exp(-i*2*pi*n*k/N);
        x2(k+1,n+1) = p;
    end
end
X = x1 * x2;
