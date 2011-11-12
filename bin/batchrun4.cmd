set nterms=300
del logfile.data
@for /L %%i in (1,1,100) do (
    submodular -m 4 -n 1000 -nterms %nterms% -fixetal -lp -heuristic -iterate
)
copy /Y logfile.data run_4_1000_%nterms%.data
