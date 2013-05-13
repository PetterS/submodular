set nterms=1
del logfile.data
@for /L %%i in (1,1,1000) do (
    submodular.exe -n 4 -m 4 -nterms  %nterms% -generators -exhaustive -generators-file generators\newgenerators.txt
)
copy /Y logfile.data run_4_1000_%nterms%.data
