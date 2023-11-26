function singleI
    for D = [5 10]
        for Run = 1:30
            path = sprintf('../IDBPI_POP/Init_Pop_Dim%02d_Run%02d.mat',D,Run);
            load(path,"-mat",'pop_I');
            pop_I = single(pop_I);
            path = sprintf('../IDBPI/Init_Pop_Dim%02d_Run%02d.mat',D,Run);
            save(path,'-mat','pop_I');
            fprintf('IDBPI DIM %d Run %d is finished.\n',D,Run);
        end
    end
end