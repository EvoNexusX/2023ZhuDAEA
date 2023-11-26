function [new_pop] = FDBPI(lb, ub, pop, val, init_popsize, dim,algRand)
    
    dim = single(dim);
    init_popsize = single(init_popsize);
    lb = single(lb);
    ub = single(ub);
    Cd  = single((pi^(dim/2))/(gamma(dim/2 + 1)));
    R   = single((((1/(init_popsize))*prod(ub - lb))/Cd)^(1/dim));
    
    h = R;
    k = 10;
    
    if min(val) < 0
        val = val - min(val);
    end
    sum_val = sum(val);

    rho = single([]);
    ppl = single([]);

    for i = single(1 : k*10)
        temp_ppl = single(lb + (ub - lb) .* rand(algRand,init_popsize/10, dim));     
        ppl = [ppl; temp_ppl];
        dist = pdist2(temp_ppl, pop);
        temp = single((1/sum_val) * sum(val' .* exp(-dist.^2 ./ (2 * h^2)), 2));
        rho = [rho; temp];  
        
        pro = i/1000;
    %     fprintf('%2d.%2d| progress: %.2f  \n', func, runs, pro);
    end
    
    [~, isort] = sort(rho, 'descend');
    isort = single(isort);
    new_pop = ppl(isort(1 : init_popsize), :);
    
    end

