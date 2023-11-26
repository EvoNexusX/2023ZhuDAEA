max_run = 10;
delete(gcp('nocreate'));
for func = 1:24
   delete(gcp('nocreate'));
   parpool('local',max_run);
   spmd(max_run)
       disp(func),disp(spmdIndex);
       Main(func, spmdIndex,10);
   end 
end
