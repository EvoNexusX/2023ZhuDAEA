function [peak,allpeak] = Main(Fn,Run,step)
  
    %% Choose the Funciton and Run
    if nargin == 0
        Fn = 4;
        Run =26; 
        step = 10;
    end
    warning('off');

    %% Showing the problem's information
    fprintf("----Function %d Run %d is running----\n",Fn,Run);
    
    %% Initialize the problem's parameters
    pro = DMMOP(Fn);
    pro1 = DMMOP1(Fn);
    peak_all = 0;

    %% Setup random seed parameters
    algRand = RandStream.create('mt19937ar','seed', Run);
    RandStream.setGlobalStream(algRand);
      
    %% Initialize the related parameters
    D = pro.D;  % The dim of the problem
    min_popsize = 7 + floor(3*log(D)); % The minimum size of the population
    lambda = min_popsize; % The minimum size of the population
    data = []; % Store the fitness
    predict_fits = single(zeros(D,step));
    correct = 0;
   
    %% The set of the global optimal solutions
    bestmem_set = [];       
    bestval_set = [];
    temp_best_pop  = [];

    %% The set of the global optimal solutions 
    old_bestmem_set = [];

    %% Setup the network
    input_size = 60-2*step;
    output_size = step;
    hidden_size = 300; 
    layers = [ ...
    sequenceInputLayer(input_size)
    lstmLayer(hidden_size)
    fullyConnectedLayer(output_size)
    regressionLayer];
    
   
    %% Test Area ---- Skip the previous environments and test the specified environment
%     while ~pro.Terminate()
%         if pro.env==40
%             break 
%         end
%         fprintf("%d\n",pro.env+1);
%         while ~pro.CheckChange(bestmem_set,bestval_set)
%             rest = pro.freq - rem(pro.evaluated, pro.freq);
%             useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
%             useless_fit = pro.GetFits(useless_pop);         
% %             if pro.CheckChange(bestmem_set,bestval_set)
% %                 continue;
% %             end
%         end
%     end 
    
    %% Initialization of the populations - Use IDBPI to generate populations based on density (IDBPI)
    path = sprintf('../IDBPI/Init_Pop_Dim%02d_Run%02d.mat',D,Run);
    load(path,"-mat",'pop_I');
    
%% Determine if the problem is terminated 
    while ~pro.Terminate()
        
        %% Get the fitness of the populations
        [fits,layers,predict_fits] = GetFitness(pro,pop_I,layers,predict_fits,Fn,Run,step);

        %% FDBPI and evaluate the fitness
        pop_F = FDBPI(pro.lower, pro.upper, pop_I, fits, 0.1*pro.freq, D,algRand);
       
        fits_F= pro.GetFits(pop_F); 
        init_pop = [pop_I;pop_F];
        fits = [fits;fits_F];
        clear fits_F pop_F;

        %% NBC_LP
        [fits, sort_index] = sort(fits, 'descend');
        init_pop = init_pop(sort_index, :);
        species = NBC(init_pop);
        species_arr = [species.len];
        imp_spec_count = length(species_arr(species_arr>=min_popsize));% The number of the important species
        [~,sort_index] = sort(species_arr,'descend');

        %% Generate groups based on clusters and record the groups's postions
        groups = init_groups(pro,lambda,init_pop,species,fits,sort_index,algRand,bestmem_set);
        clear species init_pop species species_arr fits sort_index;
        track  = track_record(groups);

        %% Initialize the archive solutions
        bestmem_set = [];
        bestval_set = [];

        %% Sort the groups by fitness
        [~,index] = sort([groups.bestval],'descend');
        groups = groups(index);
        num_groups = length(groups);
        temp_pop = zeros(num_groups*lambda,pro.D);
        % Store all the population
        OPTS = cat(1,groups.OPTS);
        temp_pop = cat(1,OPTS.pop);
        
        %% Initialize the contribution of each groups
        itermax = min(20,ceil((0.25*pro.freq)/(num_groups*min_popsize)));
        
        % Optimize all the sub-populations
        ind = 1;
        while ind<=num_groups
            [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,track);
            if ~isempty(new_groups)
                groups(ind) = new_groups;
                temp_pop((ind-1)*lambda+1:ind*lambda,:) = groups(ind).OPTS.pop;
            else
                groups(ind)= [];
                ind = ind - 1;
                num_groups = num_groups-1;
            end
            ind = ind + 1;
        end

        num_groups = length(groups);
        val = [groups.bestval];  pop = cat(1, groups.bestmem);
        [bestval, ibest] = max(val);  bestmem = pop(ibest, :);
        bestval = round(bestval);
        [~,index] = sort([groups.bestval],'descend');
        groups = groups(index);% Sort the groups by the bestmem
        OPTS = cat(1,groups.OPTS);
        temp_pop = cat(1,OPTS.pop);
        track = [track;track_record(groups)];

        %% Iteration
        while ~pro.CheckChange(double(bestmem_set),double(bestval_set))
            %% Choose the best and second-best groups
            num_groups = length(groups);
            i = 1;
            while i<=num_groups && ~isempty(groups)
                if ( groups(i).OPTS.sigma<0.02)&& ceil(groups(i).bestval)~=ceil(bestval)
                    disp(groups(i).bestval);
                    groups(i) = [];
                    i = i - 1;
                    num_groups = num_groups - 1;
                end
                i = i + 1;
            end

            % Restart if the group is empty
            if isempty(groups)
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                if rest == 0 || rest == pro.freq
                    continue;
                else
                    disp("----Restart!----");
                    [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                    continue;                
                end
            end

            fprintf("Evaluated:%d\n",rem(pro.evaluated,pro.freq));
            val = [groups.bestval]; pop = cat(1, groups.bestmem);
            [~, first_idx] = max(val);

            delta = [groups.delta];
            expected_gen = ceil((bestval - val)./(delta./itermax));
            expected_gen(first_idx) = Inf;
            if ~isempty(bestmem_set)
                gdis = pdist2(pop, bestmem_set);
                gdis = min(gdis, [], 2);
                temp_arr = [expected_gen', -gdis];
                [~, idx] = sortrows(temp_arr);
            else
                randnum = randperm(length(groups));
                temp_arr = [expected_gen', randnum'];
                [~, idx] = sortrows(temp_arr);
            end
            second_idx = idx(1);
         
            %% The best evolved populations
            ind = first_idx;
            itermax = 20;
            rest = pro.freq - rem(pro.evaluated, pro.freq);
            if pro.change == 1
                continue;
            end
            if min_popsize*itermax <= rest
                [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,track);
                track = [track;track_record(new_groups)];
            else
                itermax = min(itermax, floor(rest/min_popsize));
                if itermax>0
                    [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,track);
                    track = [track;track_record(new_groups)];
                end            
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                if rest == 0 || rest == pro.freq
                    continue;
                else
                    disp("----Restart!----");
                    [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                    track = [track;track_record(groups)];
                    continue; 
                end
            end
            
            %% Update the groups
            idx = first_idx;
            if ~isempty(new_groups)
                groups(ind) = new_groups;
                if groups(ind).bestval > bestval
                    bestmem = groups(ind).bestmem;
                    bestval = groups(ind).bestval;
                end
            else
                groups(ind)= [];
                continue;
            end
            num_groups = length(groups);
            [~,index] = sort([groups.bestval], 'descend');
            groups = groups(index);

            % Restart if the group is empty
            if isempty(groups)
                disp("----Restart!----");
                [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                track = [track;track_record(groups)];
                continue;
            end
            OPTS = cat(1,groups.OPTS);
            temp_pop = cat(1,OPTS.pop);

            %% The second-best evolved populations
            ind = second_idx;
            itermax = 20;
            rest = pro.freq - rem(pro.evaluated, pro.freq);
            if pro.change == 1
                continue;
            end
            if min_popsize*itermax <= rest
                [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,track);
                track = [track;track_record(new_groups)];
            else
                itermax = min(itermax, floor(rest/min_popsize));
                if itermax>0
                    [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,track);
                    track = [track;track_record(new_groups)];
                end            
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                if rest == 0 || rest == pro.freq
                    continue;
                else
                     
                    disp("----Restart!----");
                    [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                    track = [track;track_record(groups)];

                    continue;
                end
            end
            
            %% Update the groups
            if ~isempty(new_groups)
                groups(ind) = new_groups;
                if groups(ind).bestval > bestval
                    bestmem = groups(ind).bestmem;
                    bestval = groups(ind).bestval;
                end
            else
                groups(ind)= [];
            end
            num_groups = length(groups);
            [~,index] = sort([groups.bestval],'descend');
            groups = groups(index);
 
            % Restart if the group is empty
            if isempty(groups)
                disp("----Restart!----");
                [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                track = [track;track_record(groups)];
                continue;
            end
            OPTS = cat(1,groups.OPTS);
            temp_pop = cat(1,OPTS.pop);

            %% Save the archive solution if it converges to the global optimum solution
            if groups(idx).cc<=1e-6
                p = 1;  
                if ~isempty(bestmem_set)
                    dis_arr = pdist2(groups(idx).bestmem, bestmem_set);
                    if abs(groups(idx).bestval-bestval)<=1e-5
                        if min(dis_arr) >= 1e-3
                            bestmem_set = [bestmem_set; groups(idx).bestmem];
                            bestval_set = [bestval_set; groups(idx).bestval];
                        end
                    end
                else
                    if abs(groups(idx).bestval-bestval)<=1e-5
                        bestmem_set = [bestmem_set; groups(idx).bestmem];
                        bestval_set = [bestval_set; groups(idx).bestval];
                    end
                end
                temp_best_pop = [temp_best_pop;groups(idx).OPTS.pop];
                groups(idx) = [];
                % Restart if the group is empty
                if isempty(groups)
                    disp("----Restart!----");
                    [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                    track = [track;track_record(groups)];
                    continue;
                end
                continue
            end

            % Restart if the group is empty
            if isempty(groups)
                disp("----Restart!----");
                [groups,bestmem_set,bestval_set] = restart(track,pro,algRand,bestmem_set,bestval_set);
                track = [track;track_record(groups)];
                continue;
            end
        end
         
        peak = length(bestval_set(abs(bestval_set-max(bestval_set))<1e-5));
        peak_all = peak_all + peak;

        fprintf("---Function%d Run%d Env%d/60 Find %d个峰---\n",Fn,Run,pro.env,peak);
        disp(peak_all);
        temp_best_pop = [];
        clear track;
    end
    [peak, allpeak] = pro.GetPeak();
    PR = sum(peak, 2) ./ sum(allpeak, 2);
    disp(sum(peak,2));
    fprintf("1e-3:")
    disp(peak(1,:));
    fprintf("1e-4:")
    disp(peak(2,:));
    fprintf("1e-5:")
    disp(peak(3,:));
    disp(allpeak);
    disp(PR);

    filename = sprintf('../PEAKS/F%d_Run%d.mat', Fn, Run);
    save(filename,'peak','allpeak',"PR");

end