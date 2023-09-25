addpath(genpath(['.' filesep]));
%% setting
sig_prop = 1;
missing_p = 0;
rng(1818)
r=2;
t=10;
n=20000;
p=20;
q=2;
cor_t = eye(t);
cor_t(cor_t==0)=0.3;
cor_p = randcorr(p);
cor_sg = kron(cor_t, cor_p);

tmp = mvnrnd(0*ones(t*p,1),cor_sg,n);
Xmat = [];
for i = 1:t
    tmp1 = tmp(:, (p*(i-1)+1):(p*i));
    Xmat = [Xmat; tmp1];
end
Xmat = Xmat - mean(mean(Xmat));
X = ten_mat(Xmat, t, 't');
%
cb =0.04;
U = binornd(2,0.5,p,r);
U = mvnrnd(-2*ones(r,1), eye(r), p).*(U==1) + mvnrnd(1*ones(r,1), eye(r), p).*(U==0)...
    +mvnrnd(3*ones(r,1), eye(r), p).*(U==2);
V = mvnrnd(3*ones(r,1), 2*eye(r), t);

Z = 7*rand(n, q);
A = [1; 2];
Z = bsxfun(@minus, Z, mean(Z));


B=cb*U*V';
BZ=randsample(p, round(p*(1-sig_prop)));
B(BZ,:) = 0;

E=3*randn(n,1);
E=E-mean(E);
Y=double(ttt(tensor(X),tensor(B),[2,3],[1,2]))+Z*A+E;

var(double(ttt(tensor(X),tensor(B),[2,3],[1,2])))/var(Y)
var(Z*A)/var(Y)

%% estimate
rng(1818)
[Best,Uest,Vest,Aest,BIC]=TRIO(X, Y, 2, Z,...
    'lam1', 0, 'lam2', 0, 'fig', false, 'Niter', 5000, 'Tol', 0.01);
[Best,Uest,Vest,Aest,BIC]=TRIO(X, Y, 2, Z,...
    'lam1', 0, 'lam2', 0, 'fig', true, 'Niter', 5000, 'Tol', 0.01, 'initial', 'rand');
norm(B-Best,'fro')
norm(B,'fro')
A
Aest


BICt=[];
for lam1=linspace(0,10,10)
[Best,Uest,Vest,Aest,BIC]=TRIO(X, Y, 2, Z,...
    'lam1', lam1, 'lam2', 0, 'fig', false, 'Niter', 5000, 'Tol', 0.00001);
BICt=[BICt,BIC];
end
plot(log(BICt))

