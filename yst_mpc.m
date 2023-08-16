function P=yst_mpc(A,L)
%% 
% input:
% A- the adjacency matrix of multilayer network
% L- the number of layer

A(isnan(A)) = 0;
N = size(A,1) / L;
O = zeros(N, N);
for l = 1:L
    O = O + A(1+(l-1)*N : l*N,1+(l-1)*N : l*N);
end % layers


o = sum(O,2); % overlaping strength vector
m = L/(L-1);
P = ones(size(o)) .* m .* (o~=0); % multilayer participation coefficients

for l = 1:L

    ind = (l-1)*N+1:l*N;
    tmpA = full(A(ind,ind));
    k = sum(tmpA - diag(diag(tmpA)), 2);

    P(o~=0) = P(o~=0) ...
        - m * (k(o~=0)./o(o~=0)).^2;

end % layers