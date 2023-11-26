function [fits_I,layers,predict_fits,correct] = GetFitness(pro,pop_I,layers,predict_fits,Fun,Run,step,correct)

    %% Evaluate or predict the fitness of the populations
    if pro.env+1<60-step
        fits_I = pro.GetFits(pop_I);
        path = sprintf('../data/FUN%dRUN%d.mat',Fun,Run);
        if exist(path,'file') == 2
            load(path,'data');
            if size(data,2)==pro.env
                data = [data,fits_I];
            else
                data = fits_I;
            end
        else
            data = fits_I;
        end
        save(path,'data');
        clear data;

    elseif pro.env+1 == 60-step
        fits_I = pro.GetFits(pop_I);
        path = sprintf('../data/FUN%dRUN%d.mat',Fun,Run);
        load(path,'data');
        data = [data,fits_I];
        save(path,'data');
        mu = mean(data,2);
        sigma = std(data,0,2);
        normal_data = (data-mu)./sigma;
        options = trainingOptions('adam', 'MaxEpochs', 150, 'MiniBatchSize', 50);
        path = sprintf('../Net/FUN%dRUN%d.mat',Fun,Run);
        if exist(path) == 2
            load(path);
        else
            net = trainNetwork(normal_data(:,1:60-2*step)', normal_data(:,60-2*step+1:60-step)', layers, options);
            save(path,'net');
        end
        layers = net.Layers;
        predict_fits = predict(net,normal_data(:,step+1:60-step)')';
        predict_fits = predict_fits.*sigma+mu;
        clear data;
    else
        fits_I = pro.GetFits(pop_I(1:2000,:));
        p_fits = predict_fits(:,pro.env+1-(60-step));
        [~,index1] =  sort(fits_I);
        [~,rank1]  =  sort(index1);
        [~,index2] =  sort(p_fits(1:2000,:));
        [~,rank2]  =  sort(index2);
        sums = length(index2)*length(index2)/2;
        path = sprintf('../Net/FUN%dRUN%d.mat',Fun,Run);
        total = sum(abs(rank1-rank2));
        ac = 1-total/sums;
        correct = [correct;ac];
        fprintf("Accuracy:%.1f%%\n",ac*100);
        if ac>=0.50
            disp("Success!\n");
            fits_I = [fits_I;p_fits(2001:end,:)];
        else
            fits_I = [fits_I;pro.GetFits(pop_I(2001:end,:))];
        end
        path = sprintf('../data/FUN%dRUN%d.mat',Fun,Run);
        load(path,'data');
        data = [data,fits_I];
    end
 
%  
end