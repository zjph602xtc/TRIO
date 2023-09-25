function Y = ten_mat(X, long, dir)
%% if dir = 't' (time)
% if long = 0, then X is from a tensor (n, p, T) to matrix (nxT, p).
% if long = t, then X is from a matrix (nxt, p) to a tensor (n, p, t).
if dir == 't'
    if long == 0
        Y = double(tenmat(X,[1 3],[2]));
    else
        Y = permute(shiftdim(tensor(X', [size(X,2) size(X,1)./long long]),3),[2 1 3]);
    end
    
%% if dir = 'm' (microbiome)
% if long = 0, then X is from a tensor (n, p, T) to matrix (nxp, T).
% if long = p, then X is from a matrix (nxp, T) to a tensor (n, p, t).
elseif dir == 'm'
    if long == 0
        Y = double(tenmat(X,[1 2],[3]));
    else
        Y = shiftdim(tensor(X', [size(X,2) size(X,1)./long long]),1);
    end
%% if dir = 'n', we get X_(1), long = p
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