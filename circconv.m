function [y] = circconv(x,h,N);
% To find the circular convolution of two sequences
% x = input sequence 1
% h = impulse sequence 2
% N = Number of points in the output sequence
N2 = length(x);
N3 = length(h);
x = [x zeros(1,N-N2)]; %Append N-N2 zeros to the input sequence 1
h = [h zeros(1,N-N3)]; %Append N-N3 zeros to the sequence 2
% circular shift of the sequence 2
m = [0:1:N-1];
M = mod(-m,N);
h = h(M+1);
for n = 1:1:N
    m = n-1;
    p = 0:1:N-1;
    q = mod(p-m,N);
    hm = h(q+1);
    H(n,:) = hm;
end
% Matrix convolution
y = x*H';
