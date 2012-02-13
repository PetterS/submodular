
rem Any errors in perstistencies and lower bounds will be logged
rem by the main program into submodular.errorlog

@for /L %%i in (1,1,1000) do (
    submodular -m 4 -n 15 -nterms 10 -fixetal -optimal -heuristic -generators -exhaustive -iterate -nolog
)

