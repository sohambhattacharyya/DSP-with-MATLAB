function y = ovrlsav(x,h,N)
% To compute the output of a system using overlap and save method
% x = input seuence
% h = impulse sequence
% N = Length of each block
if(N<length(h))
    error('N must be >=length(h)')
end
Nx = length(x);
M = length(h);
M1 = M-1;
L = N-M1;
x = [zeros(1,M-1),x,zeros(1,N-1)];
h = [h zeros(1,N-M)];
K = floor((Nx+M1-1)/(L)); % Number of blocks
Y = zeros(K+1,N);
%Dividing the sequence into two blocks
for k = 0:K
    xk = x(k*L+1:k*L+N);
    Y(k+1,:) = circconv(xk,h,N);
end
Y = Y(:,M:N)'; %Discard first M-1 blocks
y = (Y(:))'
