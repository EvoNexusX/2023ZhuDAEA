%% Initialization of the populations - Use IDBPI to generate populations
max_run = 10;
delete(gcp('nocreate'));
parpool('local',max_run);
spmd(max_run)
    pro = DMMOP(1);% D is 5
    D = pro.D;
    IDBPI(pro.lower, pro.upper, 0.4*pro.freq, D, labindex);
    fprintf("Run%d is init over!",labindex);
end 
spmd(max_run)
    pro = DMMOP(1);% D is 5
    D = pro.D;
    IDBPI(pro.lower, pro.upper, 0.4*pro.freq, D, labindex+10);
    fprintf("Run%d is init over!",labindex);
end 
spmd(max_run)
    pro = DMMOP(1);% D is 5
    D = pro.D;
    IDBPI(pro.lower, pro.upper, 0.4*pro.freq, D, labindex+20);
    fprintf("Run%d is init over!",labindex);
end 


