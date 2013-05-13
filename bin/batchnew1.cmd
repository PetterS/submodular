set nterms=100
del logfile.data
@for /L %%i in (1,1,100) do (
    submodular.exe -n 1000 -m 4 -nterms  %nterms% -generators -fixetal -optimal -heuristic -iterate -m-reduction
)
copy /Y logfile.data run_4_1000_%nterms%.data
