function [B,U,V,A,varargout]=TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r, varargin)
% This is an internal function to conduct TRIO with variance/covariance
% matrix.
%    [B, U, V, A] = TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r)
%
%    [B, U, V, A] = TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r, SigmaZ, SigmaXZ, SigmaZY)
%
%    [B, U, V, A] = TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r, 'lam1', 0, 'lam2', 0, ...)
%
%    [B, U, V, A, BIC] = TRIO_corr(SigmaXY, SigmaX1, SigmaX2, r, 'n', n, 'Y', Y)
%
% -SigmaXY: C matrix, pxt, the cross-covariance matrix between the mode-1 matricization of X and y
% -SigmaX1: Theta matrix, pxp, the covariance matrix in the feature domain
% -SigmaX2: Gamma matrix, txt, the correlation matrix in the time domain
% -r: the rank of coefficent Beta
% -SigmaZ: S_zz matrix, qxq, the covariance matrix of Z
% -SigmaXZ: S_xz matrix, ptxq, cross-covariance matrices between the mode-1 matricization of X and Z
% -SigmaZY: S_zy matrix, qx1, cross-covariance matrices between Z and y
% -'n': (only necessary when request BIC) the number of sample
% -'Y': (only necessary when request BIC) raw data of Y
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

% Tianchen Xu
addpath(genpath(['.' filesep]));
[p,t]=size(SigmaXY);
inp=inputParser;
addRequired(inp,'SigmaX1',@(x)validateattributes(x,{'numeric'},{'square','ncols',p}));
addRequired(inp,'SigmaX2',@(x)validateattributes(x,{'numeric'},{'square','ncols',t}));
addRequired(inp,'r',@(x)validateattributes(x,{'numeric'},{'scalar','>=',1}));
addOptional(inp,'SigmaZ',[],@(x)validateattributes(x,{'numeric'},{'square'}));
addOptional(inp,'SigmaXZ',[]);
addOptional(inp,'SigmaZY',[]);
addParameter(inp,'lam1',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'lam2',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'lam1_ratio',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'lam2_ratio',0,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'initial','OLS',@(x)any(validatestring(x,{'OLS','rand'})));
addParameter(inp,'constraint',false,@islogical);
addParameter(inp,'rho',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
% addParameter(inp,'k',10,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
% addParameter(inp,'maxrho',inf);
addParameter(inp,'Tol',0.01,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'Niter',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0}));
addParameter(inp,'fig',false,@islogical);
if nargout > 4
    addParameter(inp,'n',[], @(x)validateattributes(x,{'numeric'},{'scalar','>=',1}));
    addParameter(inp,'Y',[], @(x)validateattributes(x,{'numeric'},{'ncols',1}));
end

parse(inp ,SigmaX1, SigmaX2, r, varargin{:});

hasz=false;
if ~isempty(inp.Results.SigmaZ)
    SigmaZ=inp.Results.SigmaZ;
    q=size(SigmaZ,1);
    SigmaXZ=inp.Results.SigmaXZ;
    SigmaZY=inp.Results.SigmaZY;
    validateattributes(SigmaXZ,{'numeric'},{'nrows',p*t,'ncols',q},'input SigmaXZ');
    validateattributes(SigmaZY,{'numeric'},{'nrows',q,'ncols',1},'input SigmaZY');
    hasz = true;
end

if nargout > 4
    n=inp.Results.n;
    n=max([n,1]);
    Y=inp.Results.Y;
    validateattributes(Y,{'numeric'},{'nrows',n},'input Y');
else
    n=1;
    Y=0;
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

% key quantities
Omega=GetOmega(t);

% initial value
switch inp.Results.initial
    case 'OLS'
        [Ub, Sb, Db]=svd(SigmaX1);
        Sb(abs(Sb)<1e-5)=0;
        B=pinv(Ub*Sb*Db)*SigmaXY*pinv(SigmaX2);
    case 'rand'
        B=rand(p,t);
end

[U,Dd,V]=svds(B,r);
U=U*Dd;
Lambda=zeros(p,t); % lagrange multiplier
if hasz
    A= pinv(SigmaZ)*SigmaZY;
    obj_ini=ObjVal_z(B,U,V,A,SigmaXY,SigmaX1,SigmaX2,SigmaZ, SigmaZY, SigmaXZ, Omega,lam1,lam2,n,Y);
else
    A = 0;
    SigmaXZ = 0;
    obj_ini=ObjVal(B,U,V,SigmaXY,SigmaX1,SigmaX2,Omega,lam1,lam2,n,Y);
end

% tsigma = zeros(inp.Results.Niter, 1);
% tsigma(1) = rho;
% for i=2:inp.Results.Niter
%     tsigma(i) = tsigma(i-1)./(sqrt(1+(1./rho).*tsigma(i-1)));
% end


%% ADMM
niter=0;
diff=inf;
rec_obj=obj_ini; % record obj value
rec_primal=[]; % record total primal residual
rec_dual=[]; % record total dual residual
rec_A=[]; % record A
while niter < inp.Results.Niter && abs(diff)>inp.Results.Tol
    niter=niter+1;
    B_old=B;
    A_old=A;

    %%%%%%%%%%%%% Primal Steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimate B
    tmp=2*SigmaXY+rho*U*V'-Lambda;
    E=tmp(:)-2*SigmaXZ*A;
    if inp.Results.constraint
        D=kron(eye(t), ones(1, p))/(2*kron(SigmaX2, SigmaX1)+rho.*eye(t*p));
        Phi = linsolve(D*kron(eye(t), ones(p,1)), D*E)';
        Phi = ones(p,1) * Phi;
    else
        Phi = 0;
    end
    Bb=(2*kron(SigmaX2, SigmaX1)+rho.*eye(t*p))\(E-Phi(:));
    B = reshape(Bb, p, t);

    % estimate U
    if lam1==0 % procrustes problem
        U=(B+(1/rho)*Lambda)*V/(V'*V);
    else % lasso
        tmp_V = mat2cell(V, t, r);
        tmp_V = repmat(tmp_V, p, 1);
        tmp_V = blkdiag(tmp_V{:});
        tmp_resp = (B + Lambda./rho)';
        tmp_resp = tmp_resp(:);
        %   [tmp_U, ~] = group_lasso(tmp_V, tmp_resp, lam1./rho, repmat(r, p, 1), 1.0, 1.0);
        opts.tFlag=5;       % run .maxIter iterations
        opts.maxIter=500;   % maximum number of iterations
        opts.nFlag=0;       % without normalization
        opts.rFlag=0;       % the input parameter 'rho'
        opts.ind=[0 r:r:(p*r)];       % set the group indices
        [tmp_U, ~, ~]= glLeastR(tmp_V, tmp_resp, lam1./rho, opts);
        U = reshape(tmp_U, r, p)';
    end

    % estimate V
    Vb = lyap((2*lam2/rho)*Omega,U'*U,(-B'-Lambda'/rho)*U);
    V = reshape(Vb, t, r);

    % normalization
    [U,D,V]=svds(U*V',r);
    U=U*D;

    % update A
    if hasz
        A = SigmaZ\(SigmaZY-SigmaXZ'*B(:));
    end

    %%%%%%%%%%%%% Dual Steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Lambda=Lambda+rho*(B-U*V');

    %%%%%%%%%%%%% Stopping Rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % primal and dual residuals
    primal=norm(B-U*V','fro')/norm(B,'fro'); % percentage difference
    rec_primal=[rec_primal,primal];
    dual=norm(B-B_old,'fro')/norm(B,'fro'); % percentage change
    rec_dual=[rec_dual,dual];
    AAAA=norm(A-A_old,'fro')/norm(A,'fro'); % percentage change
    rec_A=[rec_A,AAAA];

    % objective function value
    if hasz
        obj=ObjVal_z(B,U,V,A, SigmaXY, SigmaX1,SigmaX2,SigmaZ, SigmaZY, SigmaXZ, Omega,lam1,lam2,n,Y);
    else
        obj=ObjVal(B,U,V,SigmaXY, SigmaX1,SigmaX2,Omega,lam1,lam2,n,Y);
    end
    rec_obj=[rec_obj,obj];

    % stopping rule
    diff=(rec_obj(1,end-1)-rec_obj(1,end))./rec_obj(1,end);
    % update rho
    %     if inp.Results.varyrho
    %         rho = tsigma(max(1, 1+floor(niter/inp.Results.k)));
    %     end

    % Check Figures
    if inp.Results.fig
        % obj fcn values
        figure(101);clf;
        plot(0:niter,rec_obj,'bo-');
        title(['Objective function value (decrease in full=',num2str(rec_obj(end-1)-rec_obj(end)),')']);
        drawnow;
        % primal and dual residuals
        figure(102);clf;
        subplot(1,3,1)
        plot(1:niter,rec_primal,'o-');
        title(['|D-UV|^2: ',num2str(primal)]);
        subplot(1,3,2)
        plot(1:niter,rec_dual,'o-');
        title(['Dual residual |B-B|^2: ',num2str(dual)]);
        subplot(1,3,3)
        plot(1:niter,rec_A,'o-');
        title(['Dual residual |A-A|^2: ',num2str(AAAA)]);
        drawnow
    end
end


if niter==inp.Results.Niter
    disp(['FACTR does NOT converge after ',num2str(inp.Results.Niter),' iterations!']);
else
    disp(['FACTR converges after ',num2str(niter),' iterations.']);
end

if nargout > 4
    trueobj=obj-lam1*sum(sqrt(sum(U.^2,2)))-lam2*trace(V'*Omega*V);
    BIC = ((p+t-r)*r+sum(sum(abs(U)<1e-8))).*log(n) ...
        - 2.*(-0.5.*n.*log(2*pi)-0.5.*n.*log(trueobj./n)-n./2);
    varargout{1} = BIC;
end

end


function obj=ObjVal_z(B,U,V,A,SigmaXY,SigmaX1,SigmaX2,SigmaZ,SigmaZY,SigmaXZ,Omega,lam1,lam2,n,Y)
obj=n*(Y'*Y./n-2*trace(B'*SigmaXY)+trace(B'*SigmaX1*B*SigmaX2)...
    +2*A'*SigmaXZ'*B(:)-2*SigmaZY'*A+A'*SigmaZ*A)+lam1*sum(sqrt(sum(U.^2,2)))+lam2*trace(V'*Omega*V);
end
function obj=ObjVal(B,U,V,SigmaXY,SigmaX1,SigmaX2,Omega,lam1,lam2,n,Y)
obj=n*(Y'*Y./n-2*trace(B'*SigmaXY)+trace(B'*SigmaX1*B*SigmaX2))...
    +lam1*sum(sqrt(sum(U.^2,2)))+lam2*trace(V'*Omega*V);
end


function [Omega,VOmega,DOmega]=GetOmega(p)
% get p*p smoothing matrix Omega for smoothing spline (rank p-2)
% VOmega is p*(p-2) left singular matrix
% DOmega is (p-2)*1 singular values
% Note:
% This function only returns Omega for evenly spaced p-dimensional data
Qvv=eye(p)*(-2);
Qvv=spdiags(ones(p,1),1,Qvv);
Qvv=spdiags(ones(p,1),-1,Qvv);
Qvv=Qvv(:,2:(end-1));
Rvv=eye(p-2)*(2/3);
Rvv=spdiags(ones(p-2,1)*(1/6),1,Rvv);
Rvv=spdiags(ones(p-2,1)*(1/6),-1,Rvv);
tempv=chol(inv(Rvv))*Qvv'; % use cholesky decomp to fight against matlab precision of symmetry
Omega=tempv'*tempv; % p*p  singular matrix (rank=p-2)
[VOmega,DOmega,~]=svds(full(Omega),p-2);
DOmega=diag(DOmega);
end




