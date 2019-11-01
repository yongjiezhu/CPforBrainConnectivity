function Uinit = cp_initfactor(Y,R,init,dimorder,orthoforce)
N = ndims(Y);
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end);
        if ~isequal(size(Uinit{n}),[size(Y,n) R])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    if strcmp(init(1:6),'random')
        Uinit = cell(N,1);
        for n = dimorder
            Uinit{n} = rand(size(Y,n),R);
            Uinit{n} = bsxfun(@rdivide,Uinit{n},sqrt(sum(Uinit{n}.^2)));
        end
    elseif strcmp(init(1:5),'nvecs') || strcmp(init(1:4),'eigs')
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            fprintf('Computing %d leading e-vectors for factor %d.\n',R,n);
            Uinit{n} = nvecs(Y,n,R);
        end
    else
        error('The selected initialization method is not supported');
    end
end

lambda = ones(R,1);
lamembind = N;
switch orthoforce
    case 1 % 3-D tensor
        if N ==3
            [Uinit{1},Uinit{2},Uinit{3}]=dtld(double(Y),R);
            lambda_n = zeros(N-1,R);
            ind = 0;
            for n = dimorder([1:lamembind-1 lamembind+1:end])
                ind = ind +1;
                lambda_n(ind,:) = sqrt(sum(Uinit{n}.^2));
                Uinit{n} = bsxfun(@rdivide,Uinit{n},sqrt(sum(Uinit{n}.^2)));
            end
            Uinit{lamembind} = bsxfun(@times,Uinit{lamembind},prod(lambda_n));
        end
        
    case 2 %ALS init
        for n = dimorder(1:end)
            % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
            Unew = mttkrp(Y,Uinit,n);
            
            % Compute the matrix of coefficients for linear system
            W = ones(R,R);
            for i = [1:n-1,n+1:N]
                W = W .* (Uinit{i}'*Uinit{i});
            end
            Unew = (W \ Unew')';
            Unew = max(eps,(Unew));
            % Normalize each vector to prevent singularities in coefmatrix
            lambda = sqrt(sum(Unew.^2,1))'; %2-norm
            Unew = Unew * spdiags(1./lambda,0,R,R);
            Uinit{n} = Unew;
        end
        Uinit{end} = bsxfun(@times,Uinit{end},lambda');
end