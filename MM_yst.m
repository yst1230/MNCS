%%%graph learning via MM

%%Inputs:
%X_noisy: noisy data matrix of size p*nSamples
%a and b: the hyperparameters alpha and beta
%w_0: initial value of weight vector w
%thresh: Convergence threshold to exit the loop
%nSamples: Number of realizations of graph signal

%Outputs:
%wk: the estimated weight vector
%stat: A structure array containing percentage of non-zero entities in %w,time, number of iterations and f(w) as fields.

function [wk, stat] = MM_yst(data, a, b,w_0,thresh)
%%  data weighting
[pointNum,nodeNum]=size(data);
V0=(1/pointNum)*eye(pointNum);
tmp=pointNum*(V0*data);
tmp=tmp-repmat(mean(tmp),pointNum,1); % centrlization
%% pairwise distance
D = sparse(gsp_distanz(tmp).^2);    % pairwise distance matrix
D=D./max(D(:));
d = squareform_sp(D);                 % vectored Z
m = length(d);                        % m=(p*(p-1))/2;
[S, ~] = sum_squareform(nodeNum);    % S is a binary matrix such that Sw=W1, where W is the weight matrix; St is the transpose of S
wk=w_0;
wk(wk==0)=eps;
obj_val(1) = (2*wk'*d) - (a*sum(log(S*wk)))  + (b*(norm(wk))^2) ; % objective function value
% n_z(1)=length(find(wk~=0))*100/m;                                 % Percentage of non-zero (active) elements in the weight vector w
idx=2;
eta=1;
d2=-2*d;
d4=4*(d.^2);
b8=8*b;
b4=4*b;
Swk=S*wk;
Swki=1./Swk;
Swkii=repmat(Swki,[1 m]);
SS=zeros(nodeNum,m);
%% fix v and estimate w
while eta>thresh
    SS(S~=0)=Swkii(S~=0);
    cl=a*sum(SS)'.*wk;
    u=((d2)+sqrt((d4)+(b8*cl)))/(b4);
    wk=u;
    wk(u<0.001)=0;
%     n_z(idx)=length(find(wk~=0))*100/m;
    Swk=S*wk;
    Swki=1./Swk;
    Swkii=repmat(Swki,[1 m]);
    obj_val(idx) = (2*u'*d) - (a*sum(log(S*u))) + (b*(norm(u))^2);
    eta=abs(obj_val(idx-1)-obj_val(idx));
    idx=idx+1;
end

stat.num_itr = length(obj_val); %total number of iterations
stat.obj_val = obj_val; %objective function value
end

