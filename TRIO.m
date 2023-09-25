function [B,U,V,A,BIC]=TRIO(X,Y,r,varargin)
% This is a function to conduct TRIO regression
%    [B, U, V, A, BIC] = TRIO(X, Y, r)
%
%    [B, U, V, A, BIC] = TRIO(X, Y, r, Z)
%
%    [B, U, V, A, BIC] = TRIO(X, Y, r, Z, 'lam1', 0, 'lam2', 0, ...)
%
%
% -X: centered predictor, nxpxT
% -Y: centered response, nx1
% -r: the rank of coefficent Beta
% -Z: centered adjusting covariate, nxq
% -'initial' ('OLS'): choose the initial value of coefficient Beta
%      'OLS': OLS method, B=inv(SigmaX1)*SigmaXY*pinv(SigmaX2)
%      'rand': a random matrix
% -'lam1' (0): lambda1
% -'lam2' (0): lambda2
% -'rho' (100): (only used when lam1_ratio~=0 or lam2_ratio~=0) a common parameter to determine lamdba1 and lambda2:
% lambda1=lam1_ratio*rho, lambda2=lam2_ratio*rho
% -'lam1_ratio' (0): lambda1=lam1_ratio*rho
% -'lam2_ratio' (0): lambda2=lam2_ratio*rho
% -'constraint' (false): whether include the compositional constraint (1^TB=0)
% -'Tol': percentage of likelihood change to terminate the function
% -'Niter': max iterations
% -'fig': whether show diagnosis figures

% Tianchen Xu, Gen Li, Kun Chen
addpath(genpath(['.' filesep]));

if (isnumeric(X))
    X = tensor(X);
end

X = ten_mat(X, 0, 't');

[n,~]=size(Y);
[tmp, p]=size(X);
t = tmp/n;

inp=inputParser;
addRequired(inp,'r',@(x)validateattributes(x,{'numeric'},{'scalar','>=',1}));
addOptional(inp,'Z',[],@(x)validateattributes(x,{'numeric'},{'nrows',n}));
addParameter(inp,'lam1',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'lam2',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'lam1_ratio',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'lam2_ratio',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'rho',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'initial','OLS',@(x)any(validatestring(x,{'OLS','rand'})));
addParameter(inp,'constraint',false,@islogical);
addParameter(inp,'Tol',0.01,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'Niter',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'fig',false,@islogical);

parse(inp, r, varargin{:});

Xten = ten_mat(X, t, 't');
Xmat = ten_mat(Xten, 0, 'n');

hasz=false;
if ~isempty(inp.Results.Z)
    Z=inp.Results.Z;
    [~, q]=size(Z);
    SigmaXZ=cov([Xmat Z],'partialrows');
    SigmaXZ=SigmaXZ(1:(p*t),(end-q+1):end);
    SigmaZY=cov([Z Y], 'partialrows');
    SigmaZY=SigmaZY(1:(end-1),end);
    SigmaZ=cov(Z, 'partialrows');
    hasz = true;
end

Sigma = cov(Xmat, 'partialrows');
[SigmaX2,SigmaX1,~] = NearestKroneckerProduct(Sigma,[t t],[p,p],true);
% kron.em function for Expectation-Maximization (EM) algorithm to estimate
% the separable covariance can be obtained via
% https://www.sciencedirect.com/science/article/pii/S0047259X17300477
% (appendix D)
% https://doi.org/10.1016/j.jmva.2018.03.010

SigmaXY=cov([Xmat Y], 'partialrows');
SigmaXY=SigmaXY(1:(end-1),end);
SigmaXY=reshape(SigmaXY, p, t);

if (max(abs(mean(Xmat)))>1e-5)
    warning('X not centered!')
end
if (abs(mean(Y))>1e-5)
    warning('Y not centered!')
end
if hasz && (max(abs(mean(Z)))>1e-5)
    warning('Z not centered!')
end

rho=inp.Results.rho;
if inp.Results.lam1_ratio > 0
    lam1=inp.Results.lam1_ratio.*rho;
else
    lam1=inp.Results.lam1;
end

if inp.Results.lam2_ratio > 0
    lam2=inp.Results.lam2_ratio.*rho;
else
    lam2=inp.Results.lam2;
end

if hasz
    [B,U,V,A,BIC]=TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r, SigmaZ, SigmaXZ, SigmaZY, ...
        'lam1', lam1, 'lam2', lam2, 'constraint', inp.Results.constraint,  ...
        'Tol', inp.Results.Tol, 'Niter', inp.Results.Niter, 'fig', inp.Results.fig, ...
        'n', n, 'Y', Y, 'initial', inp.Results.initial);
else
    [B,U,V,A,BIC]=TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r, ...
        'lam1', lam1, 'lam2', lam2, 'constraint', inp.Results.constraint,  ...
        'Tol', inp.Results.Tol, 'Niter', inp.Results.Niter, 'fig', inp.Results.fig, ...
        'n', n, 'Y', Y, 'initial', inp.Results.initial);
end

end

function [ B, C, D ] = NearestKroneckerProduct( A, SizeB, SizeC, Hermitian ) %#codegen
% https://doi.org/10.1007/978-94-015-8196-7_17
m = size( A, 1 );
n = size( A, 2 );
m1 = SizeB( 1 );
n1 = SizeB( 2 );
m2 = SizeC( 1 );
n2 = SizeC( 2 );
if nargin < 4
    Hermitian = false;
end
assert( m == m1 * m2 );
assert( n == n1 * n2 );
if Hermitian
    assert( m1 == n1 );
    assert( m2 == n2 );
    A = 0.5 * ( A + A' );
end
R = reshape( permute( reshape( A, [ m2, m1, n2, n1 ] ), [ 2 4 1 3 ] ), m1 * n1, m2 * n2 );
[ B, S, C ] = svds( R, 1 );
SqrtS = sqrt( S );
B = reshape( B * SqrtS, m1, n1 );
C = reshape( C * SqrtS, m2, n2 );
if Hermitian
    B = 0.5 * ( B + B' );
    C = 0.5 * ( C + C' );

    if all( diag( B ) < 0 ) && all( diag( C ) < 0 )
        B = -B;
        C = -C;
    end
end
if nargout > 2
    D = A - kron( B, C );
    if Hermitian
        D = 0.5 * ( D + D' );
    end
end
end

function Y = ten_mat(X, long, dir) % change the shape of a tensor to matrix (vice versa)
% if dir = 't' (time)
% if long = 0, then X is from a tensor (n, p, T) to matrix (nxT, p).
% if long = t, then X is from a matrix (nxt, p) to a tensor (n, p, t).
if dir == 't'
    if long == 0
        Y = double(tenmat(X,[1 3],[2]));
    else
        Y = permute(shiftdim(tensor(X', [size(X,2) size(X,1)./long long]),3),[2 1 3]);
    end

    % if dir = 'm' (microbiome)
    % if long = 0, then X is from a tensor (n, p, T) to matrix (nxp, T).
    % if long = p, then X is from a matrix (nxp, T) to a tensor (n, p, t).
elseif dir == 'm'
    if long == 0
        Y = double(tenmat(X,[1 2],[3]));
    else
        Y = shiftdim(tensor(X', [size(X,2) size(X,1)./long long]),1);
    end
    % if dir = 'n', we get X_(1), long = p
    % if long = 0, then X is from a tensor (n, p, T) to matrix (n, pxT).
    % if long = p, then X is from a matrix (n, pxT) to a tensor (n, p, t).
elseif dir =='n'
    if long == 0
        Y = double(tenmat(X, [2 3],[1]))';
    else
        Y = tensor(X, [size(X,1) long size(X,2)./long]);
    end
else
    error('wrong dir')
end

end