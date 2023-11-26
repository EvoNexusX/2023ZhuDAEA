function [group,track] = KE_CMA_ES(pro, group, lb, ub, itermax, algRand,temp_pop,idx,temp_best_pop,track)
    
    xmean   = group.xmean;
    bestmem = group.bestmem;
    bestval = group.bestval;
    OPTS    = group.OPTS;

    old_bestval = group.bestval;
    lambda = 7 + floor(3*log(pro.D));

    temp_pop((idx-1)*lambda+1:end,:)=[];
    temp_pop = [temp_pop;temp_best_pop];

    % 初始化系数
    dim = length(xmean);
    sigma = OPTS.sigma;
    
    if OPTS.first == 1
        lambda = 7 + floor(3*log(dim));
%         lambda = length(group);
        mu = floor(lambda/2);
        % Strategy parameter setting: Selection
        weights = log(mu+1/2)-log(1:mu)';       % muXone recombination weights
        mu = floor(mu);                         % number of parents/points for recombination
        weights = weights/sum(weights);         % normalize recombination weights array
        mueff=sum(weights)^2/sum(weights.^2);   % variance-effective size of mu
        
        % Strategy parameter setting: Adaptation
        cc = (4+mueff/dim) / (dim+4 + 2*mueff/dim);     % time constant for cumulation for C
        cs = (mueff+2)/(dim+mueff+5);                   % t-const for cumulation for sigma control
        c1 = 2 / ((dim+1.3)^2+mueff);                   % learning rate for rank-one update of C
        cmu = 2 * (mueff-2+1/mueff) / ((dim+2)^2+2*mueff/2);    % and for rank-mu update
        damps = 1 + 2*max(0, sqrt((mueff-1)/(dim+1))-1) + cs;   % damping for sigma
        
        % Initialize dynamic (internal) strategy parameters and constants
        pc = zeros(dim,1); ps = zeros(dim,1);           % evolution paths for C and sigma
        B = eye(dim);                                   % B defines the coordinate system
        D = eye(dim);                                   % diagonal matrix D defines the scaling
        C = B*D*(B*D)';                                 % covariance matrix
        chiN=dim^0.5*(1-1/(4*dim)+1/(21*dim^2));        % expectation of ||N(0,I)|| == norm(randn(N,1))
        countval = 0;
        iters = 0;
    else
       
        lambda = OPTS.lambda;
        weights = OPTS.weights;
        mu = OPTS.mu;
        mueff = OPTS.mueff;
        cc = OPTS.cc;
        cs = OPTS.cs;
        c1 = OPTS.c1;
        cmu = OPTS.cmu;
        damps = OPTS.damps;
        pc = OPTS.pc;
        ps = OPTS.ps;
        B = OPTS.B;
        D = OPTS.D;
        C = OPTS.C;
        chiN = OPTS.chiN;
        countval = OPTS.countval;
        iters = group.iters;
    end
    
    % -------------------- Generation Loop --------------------------------
    stopiters = iters + itermax; 
    % Fes = 0;
    while iters < stopiters
        % Generate and evaluate lambda offspring

        rho_x = zeros(1,lambda);
        temp_x =zeros(dim, 10);

        temp_arz = randn(algRand, dim, lambda);
        for j = 1 : lambda
            temp_x(:,j) = xmean + sigma*(B*D*temp_arz(:,j));
            temp_ub = temp_x(:, j) > ub(1);
            temp_lb = temp_x(:, j) < lb(1);
            if any(temp_ub) || any(temp_lb)
                temp_x(temp_ub, j) =  ub(1);
                temp_x(temp_lb, j) =  lb(1);
                temp_arz(:, j) = pinv(D) * pinv(B) * ((temp_x(:, j) - xmean)/sigma);
            end
            % countval = countval + 1;
        end
        arx = temp_x;
        arz = temp_arz;


        % ND-sort
        obj = zeros(2,lambda);
        rho = kernel(arx',temp_pop);
        group.OPTS.count =group.OPTS.count + length(rho(rho>0));

        fit = pro.GetFits(arx');
        obj(1,:) = fit';
        obj(2,:) = rho;

        [index,~] = Fast_ND_SORT(obj);

        arx = arx(:,index);
        obj = obj(:,index);
        arz = arz(:,index);
        fit = fit(index,:);

        x = group.mean_distance;

        if group.iters>0 && group.OPTS.count > 100/x
            disp("----The goup is useless!----!3\");
            disp(abs(bestval));
            group = [];
            break;
        end

        arfitness = fit;

        iters = iters + 1;

        if size(arfitness,1)~=size(arx',1)
            break;
        end

        % Sort by fitness and compute weighted mean into xmean
        arindex = 1:lambda;
        xmean = arx(:,arindex(1:mu))*weights;
        zmean = arz(:,arindex(1:mu))*weights;
      
        % Cumulation: Update evolution paths
        ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B * zmean);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*countval/lambda))/chiN < 1.4+2/(dim+1);
        pc = (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * (B*D*zmean);
        
       

        % Adapt covariance matrix C
        C = (1-c1-cmu) * C ...
            + c1 * (pc*pc' ... % plus rank one update
            + (1-hsig) * cc*(2-cc) * C) ... % minor correction
            + cmu ... % plus rank mu update
            * (B*D*arz(:,arindex(1:mu))) ...
            * diag(weights) * (B*D*arz(:,arindex(1:mu)))';
        
        % Adapt step-size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        
        delta = bestval-old_bestval;
 
        % Update B and D from C
        C = triu(C) + triu(C,1)';       % enforce symmetry
        [B,D] = eig(C);                 % eigen decomposition, B==normalized eigenvectors
        D = diag(sqrt(diag(abs(D))));   % D contains standard deviations now
        % Break, if fitness satisfies stop condition
        
        if OPTS.first == 1
            OPTS.first = 0;
        end
        
        [arfitness, index] = sort(arfitness,'descend');         % maxmization

        if arfitness(1) > bestval
            bestmem = arx(:,index(1))';
            bestval = arfitness(1);
        end
        
        if std(arfitness) < 1e-6
            break;
        end
    end
    if ~isempty(group)
        OPTS.pc = pc;
        OPTS.ps = ps;
        OPTS.B = B;
        OPTS.D = D;
        OPTS.C = C;
        OPTS.sigma = sigma;
        OPTS.lambda = lambda;
        OPTS.weights = weights;
        OPTS.mu = mu;
        OPTS.mueff = mueff;
        OPTS.cc = cc;
       
        OPTS.cs = cs;
        OPTS.c1 = c1;
        OPTS.cmu = cmu;
        OPTS.damps = damps;
        OPTS.chiN = chiN;
        OPTS.countval = countval;
        OPTS.pop = arx';
        group.xmean    = xmean;
        group.bestmem  = arx(:,1)';
        group.bestval  = arfitness(1);
        group.OPTS     = OPTS;
       
        group.delta    = bestval - old_bestval;
        group.cc       = std(arfitness);
        group.iters    = iters;
    end
end