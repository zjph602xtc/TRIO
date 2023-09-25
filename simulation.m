% Simulation
% addpath(genpath(['.' filesep]));
addpath(genpath(['.' filesep]));
allres = cell(5,12);
nullres = zeros(5,12);
rrow=1; rcol=1;
for sig_prop = [0.4]
    for smooth_p = [1]
        for missing_p = [0.5]
            rng(1818)
            % generate data
            for iter = 1:100
                r=3;
                t=10;
                %                 t=30;
                n=500;
                p=100;

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
                X = ten_mat(Xmat, t, 't');

                snr = [];
                if ~exist('cb')
                    cb =0.04;
                end
                if iter==1
                    for ss=1:100000
                        U = binornd(2,0.5,p,r);
                        U = mvnrnd(-2*ones(r,1), eye(r), p).*(U==1) + mvnrnd(1*ones(r,1), eye(r), p).*(U==0)...
                            +mvnrnd(3*ones(r,1), eye(r), p).*(U==2);
                        V = mvnrnd(3*ones(r,1), 2*eye(r), t);

                        for i=1:r
                            V(:,i)=smooth(1:t, V(:,i),smooth_p,'lowess');
                        end

                        B=cb*U*V';
                        BZ=randsample(p, round(p*(1-sig_prop)));
                        BNZ=setdiff(1:p,BZ);
                        B(BZ,:) = 0;

                        E=3*randn(n,1);
                        E=E-mean(E);
                        Y=double(ttt(tensor(X),tensor(B),[2,3],[1,2]))+E;

                        snr = [snr var(double(ttt(tensor(X),tensor(B),[2,3],[1,2])))/var(Y)];
                        if length(snr)>31
                            msnr=mean(snr((length(snr)-30):end));
                            if msnr < 0.85
                                cb=cb+0.001;
                                snr = [];
                            elseif  msnr > 0.9
                                cb=max(0.001, cb - 0.001);
                                snr = [];
                            else
                                break
                            end
                        end

                        cb;
                        mean(snr);
                    end
                else
                    E=3*randn(n,1);
                    E=E-mean(E);
                    Y=double(ttt(tensor(X),tensor(B),[2,3],[1,2]))+E;
                end
                tmp = shiftdim(repmat(binornd(1, missing_p, t, n, 1), 1, 1, p), 1);
                X = double(X);
                X(logical(tmp)) = nan;
                X = tensor(X);
                Xmat = ten_mat(X, 0, 't');

                writematrix(Xmat, ['G:\tensor_data\fig2\' sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_%d.csv',sig_prop, smooth_p, missing_p, iter)])
                writematrix(Y, ['G:\tensor_data\fig2\' sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_y_%d.csv',sig_prop, smooth_p,missing_p, iter)])

                %
                if (iter==1)
                    n=5000; % test set
                    tmp = mvnrnd(0*ones(t*p,1),cor_sg,n);
                    Xmat_t = [];
                    for i = 1:t
                        tmp1 = tmp(:, (p*(i-1)+1):(p*i));
                        Xmat_t = [Xmat_t; bsxfun(@minus, tmp1, mean(tmp1,2))];
                    end
                    X_t = ten_mat(Xmat_t, t, 't');
                    Y_t=double(ttt(tensor(X_t),tensor(B),[2,3],[1,2]));

                    writematrix(Xmat_t, ['G:\tensor_data\fig2\' sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_%d.csv',sig_prop, smooth_p, missing_p, iter)])
                    writematrix(Y_t, ['G:\tensor_data\fig2\' sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_y_%d.csv',sig_prop, smooth_p,missing_p, iter)])
                end
            end


            %% estimation
            error_ols=[];
            error_cp=[];
            error_trio=[];
            X_te = csvread(['G:\tensor_data\fig2\' sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_%d.csv',sig_prop, smooth_p, missing_p, 1)]);
            X_te_t = ten_mat(X_te, t, 't');
            Y_te = csvread(['G:\tensor_data\fig2\' sprintf('t_sig_p_%.2f_Vcor_%d_missing_%.2f_y_%d.csv',sig_prop, smooth_p,missing_p, 1)]);
            lam11=[];lam22=[];
            for iter=1:100
                iter;
                Xmat = csvread(['G:\tensor_data\fig2\' sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_%d.csv',sig_prop, smooth_p, missing_p, iter)]);
                Y = csvread(['G:\tensor_data\fig2\' sprintf('sig_p_%.2f_Vcor_%d_missing_%.2f_y_%d.csv',sig_prop, smooth_p,missing_p, iter)]);
                X = ten_mat(Xmat, t, 't');

                SigmaX1 = cov(Xmat,'partialrows');
                SigmaX2 = cov(ten_mat(X, 0, 'm'), 'partialrows');

                SigmaXY = zeros(p, t);
                for i=1:p
                    for j=1:t
                        tmp = cov(double(X(:,i,j)), Y, 'partialrows');
                        SigmaXY(i,j)=tmp(1,2);
                    end
                end

                % tensor OLS (random design) 
                B_OLS=SigmaX1\SigmaXY/SigmaX2;
                error_ols = [error_ols norm(Y_te - double(ttt(tensor(X_te_t),tensor(B_OLS),[2,3],[1,2])),'fro')/sqrt(size(Y_te,1))];


                % TRIO (proposed)
                if iter<=10
                    tmp=[];
                    l1try = 0:0.3:10;
                    l2try = 0:5:30;
                    for l1 = l1try
                        tmp1 = [];
                        for l2 = l2try
                            [B_FACTR,U_FACTR,V_FACTR,A,BIC]=TRIO_corr(SigmaXY,SigmaX1,SigmaX2,r,...
                                'fig',false,'Niter',2000,'lam1',l1,'lam2',l2,'n',n,'Y',Y);
                            tmp1 = [tmp1 norm(Y_te - double(ttt(tensor(X_te_t),tensor(B_FACTR),[2,3],[1,2])),'fro')/sqrt(size(Y_te,1))];
                        end
                        tmp = [tmp; tmp1];
                    end
                    [rrrr1,rrrr2]=find(min(min(tmp))==tmp);
                    lam1 = l1try(rrrr1);
                    lam2 = l2try(rrrr2);
                    if lam1==10
                        warning('bad lam1')
                    end
                    if lam2==30
                        warning('bad lam2')
                    end
                    lam11=[lam11 lam1];
                    lam22=[lam22 lam2];
                end
                [B_FACTR,U_FACTR,V_FACTR,A,BIC]=TRIO_corr(SigmaXY,SigmaX1,SigmaX2,r,...
                     'fig',false,'Niter',2000,'lam1',mean(lam11),'lam2',mean(lam22),'n',n,'Y',Y);
                error_trio = [error_trio norm(Y_te - double(ttt(tensor(X_te_t),tensor(B_FACTR),[2,3],[1,2])),'fro')/sqrt(size(Y_te,1))];
             
                % RegTen (from TensorReg)
                % CP regression
                try
                    [~,B_CP,~,~] = kruskal_reg([],permute(X,[2,3,1]),Y,r,'normal'); % CP low-rank regression
                    error_cp = [error_cp norm(Y_te - double(ttt(tensor(X_te_t),tensor(B_CP),[2,3],[1,2])),'fro')/sqrt(size(Y_te,1))];
                end
          
                % regularized regression (robust! even for small n)
                % [~,B_rCP,~] = matrix_sparsereg([],permute(X,[2,3,1]),Y,0,'normal'); % tuning need to be selected for nuclear penalty
            end
            allres{rrow, rcol} = [error_trio; error_cp; error_ols];
            nullres(rrow,rcol) = norm(Y_te ,'fro')/sqrt(size(Y_te,1));

            if rrow == 4
                rrow=0;
                rcol=rcol+1;
            end
            rrow = rrow +1;
        end
    end
end