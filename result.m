for Fn = 1:24
    PRS = 0;
    for Run = 1:30
        path = sprintf('PEAKS/F%d_Run%d.mat',Fn,Run);
        load(path,'PR');
        PRS = PRS + PR;
    end
    PRS = PRS/30;
    disp(PRS);
    disp('----------');
end