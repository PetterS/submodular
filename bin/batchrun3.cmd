set nterms=1000

del logfile.data
@for /L %%i in (1,1,100) do submodular -m 3 -n 1000 -nterms %nterms% -iccv11 -lp -heuristic -iterate
copy /Y logfile.data run_3_1000_%nterms%.data
