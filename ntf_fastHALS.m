function [P,U,fitarr] = ntf_fastHALS(Y,R,opts)
%PARAFAC decomposition - NTF.
%
% Ref:  Multi-way Nonnegative Tensor Factorization Using Fast
%Hierarchical Alternating Least Squares Algorithm
% Anh Huy PHAN, Andrzej CICHOCKI 04/2008


eps = 1e-100;

%% Fill in optional variable
if ~exist('opts','var')
    opts = struct;
end
%% Extract number of dimensions and norm of Y.
N = ndims(Y);

%% Set algorithm parameters from input or by using defaults
defoptions = struct('tol',1e-4,'maxiters',500,...
    'init','random','orthoforce',0,'fitmax',1,...
    'verbose',0,'lsparse',zeros(1,N),'lsmooth',zeros(1,N),...
    'lsorth',zeros(1,N),'dimorder',1:N);
if ~exist('opts','var')
    opts = struct;
end
[tol,maxiters,init,orthoforce,fitmax,...
    verbose,lsparse,lsmooth,lsorth,dimorder] = scanparam(defoptions,opts);

% dimorder = 1:N;
if isfield(opts,'dimorder')
    dimorder =  opts.dimorder;
end
normY = norm(Y)^2;

%% Set up for iterations - initializing U and the fit.
U = cp_initfactor(Y,R,init,dimorder,orthoforce);

lambda = ones(R,1);
In = size(Y);

lamembind = N;
P = ktensor(U);
normresidual = ( normY^2 + norm(P)^2 - 2 * innerprod(Y,P) );
fit = 1 - (normresidual / normY); %fraction explained by model
if verbose
    fprintf('\nfast HALS NTF:\n');
    fprintf(' Iter %2d: fit = %e \n', 0, fit );
end
fitarr = zeros(1,maxiters);
fitarr(1) = fit;

%% matricize Y to mode k
for n = 1: N
    Ymode{n} = permute(Y,[n 1:n-1,n+1:N]);
    Ymode{n} = reshape(Ymode{n}.data, In(n), []);
end

T31 = U{1}'*U{1};
for n = 2: N
    T31 = T31.*(U{n}'*U{n});
end

Uf = U;
for n =1:N
    Uf{n}(:) = 0;
end

%% Main Loop: Iterate until convergence
for iter = 1:maxiters
    fitold = fit;
    lambda = (sum(U{lamembind }.*U{lamembind }))';
    % Iterate over all N modes of the tensor
    for ksm = 1:N
        Uf{ksm} = lsmooth(ksm) * [U{ksm}(2,:); (U{ksm}(1:end-2,:)+U{ksm}(3:end,:))/2; U{ksm}(end-1,:)] ;
    end
    
    for n = dimorder(1:end)
        Uorth = 0;
        T2 = Ymode{n} * khatrirao(U([1:n-1 n+1:N]),'r') ;
        T32 = T31./((U{n}'*U{n}));
        for ri = 1:R
            if lsorth(n)
                %                  Uorth = sum(U{n}(:,[1:ri-1 ri+1:R]),2);
                Uorth = U{n}*U{n}';
                Uorth = (bsxfun(@rdivide,Uorth,diag(Uorth))- eye(In(n)))*U{n}(:,ri);
            end
            
            if n ~= lamembind
                U{n}(:,ri) = U{n}(:,ri)* lambda(ri)+ T2(:,ri) - ...
                    U{n} * T32(:,ri) -lsparse(n) + Uf{n}(:,ri) + ...
                    lsorth(n)*Uorth;   % lsorth(n)*sum(U{n}(:,[1:ri-1 ri+1:R]),2));
                
                U{n}(:,ri) = max(eps,( U{n}(:,ri)));
                
                lambda_ri = norm(U{n}(:,ri)); %2-norm
                U{n}(:,ri) = U{n}(:,ri)/lambda_ri;
            else
                U{n}(:,ri) = U{n}(:,ri) + T2(:,ri) - U{n} * T32(:,ri) + ...
                    Uf{n}(:,ri) - lsparse(n)+ ...
                    lsorth(n)*Uorth; % lsorth(n)*sum(U{n}(:,[1:ri-1 ri+1:R]),2));
                U{n}(:,ri) = U{n}(:,ri)/(1+ lsmooth(N));
                
                U{n}(:,ri) = max(eps,( U{n}(:,ri)));
                
            end
        end
        
        T31 = T32 .* (U{n}'*U{n});
    end
    
    if  ((mod(iter,10) == 1) || (iter == maxiters))
        P = ktensor(U);
        normresidual = normY + norm(P)^2 - 2 * innerprod(Y,P);
        fit = 1 - (normresidual / normY); %fraction explained by model
        fitchange = abs(fitold - fit);
        fitarr(iter) = fit;
        if verbose == 1
            fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
        end
        if (iter > 1)  && ((fitchange < tol))
            break;
        end
        
    end
end

%% Clean up final result
fitarr(iter+1:end) = [];
end