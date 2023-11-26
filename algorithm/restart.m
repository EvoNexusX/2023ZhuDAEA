function [groups,bestmem_set,bestval_set] = restart(track,pro, algRand,bestmem_set,bestval_set)
    rest = pro.freq - rem(pro.evaluated, pro.freq);
    if rest == 0 || rest==pro.freq
        groups = [];
        return;
    end
    groups = [];
    lambda =  7 + floor(3*log(pro.D));
    if rest >= lambda*2
        useless_pop = single(rand(10*rest,pro.D) .* (pro.upper - pro.lower) + pro.lower);
        rho = kernel(useless_pop,track);
        [~,index] = sort(rho);
        get_pop = useless_pop(index(1),:);
        D= pro.D;
        min_popsize = 7 + floor(3*log(D));
        groups = struct();
        i = 1;
        groups(i).idx = i;
        groups(i).OPTS.first = 1;
 
        groups(i).xmean = get_pop';
        groups(i).OPTS.pop = groups(i).xmean' + 0.5*randn(algRand, min_popsize, D);
        groups(i).OPTS.val = pro.GetFits(groups(i).OPTS.pop);
        groups(i).OPTS.count = 0;
        groups(i).OPTS.sigma = 0.5;
        groups(i).cc = std(groups(i).OPTS.val);
        [~,index] = sort(groups(i).OPTS.pop,"descend");
        groups(i).bestval = groups(i).OPTS.val(index(1));
        groups(i).bestmem = groups(i).OPTS.pop(index(1),:);
        groups(i).delta = 0;
        groups(i).iters = 0;
        groups(i).mean_distance = 10;
    else
        useless_pop = single(rand(10*rest,pro.D) .* (pro.upper - pro.lower) + pro.lower);
        rho = kernel(useless_pop,track);
        [~,index] = sort(rho);
        get_pop = useless_pop(index(1:rest),:);
        get_fit = pro.GetFits(get_pop);
        for i = 1:size(get_pop,1)
 
                if abs(get_fit(i)-max(bestval_set))<=1e-3
                    bestmem_set = [bestmem_set;get_pop(i,:)];
                    bestval_set = [bestval_set;get_fit(i,:)];
                end
 
        end
    end
 
end