function y = ovrladd(x,h,L)
% To compute the output of a system using overlap and add method
% x = input sequence
% h = impulse sequence
% L = Length of each block
Nx = length(x);
M = length(h);
M1 = M-1;
R = rem(Nx,L);
N = L+M1;
x = [x zeros(1,L-R)];
h = [h zeros(1,N-M)];
K = floor(Nx/L); % Number of blocks
Y = zeros(K+1,N);
z = zeros(1,M1);
% Dividing the sequence into K blocks
for k = 0:K
    xp = x(k*L+1:k*L+L);
    xk = [xp z];
    y(k+1,:) = circconv(xk,h,N);
end
yp = y';
[x,y] = size(yp);
for i = L+1:x
    for j=1:y-1
        temp1 = i-L;
        temp2 = j+1;
        temp3 = yp(temp1,temp2)+yp(i,j);
        yp(temp1,temp2) = yp(i,j);
        yp(temp1,temp2) = temp3;
    end
end
z = 1;
for j = 1:y
    for i = 1:x
        if((i<=L && j<=y-1)||(j==y))
            ypnew(z) = yp(i,j);
            z = z+1;
        end
    end
end
y = ypnew
