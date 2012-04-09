@echo off
echo Build file for Windows 32/64-bit. Execute from a Visual Studio Command Prompt.
echo Might need to be modified to find Python include/libs directories
echo Petter Strandmark 2012

set OLDLIB=%LIB%
set OLDINCLUDE=%INCLUDE%

set INCLUDE=../thirdparty/maxflow-v3.01.src/;%INCLUDE%
set INCLUDE=../thirdparty/Petter/;%INCLUDE%
set INCLUDE=../source/library/;%INCLUDE%

set LIB=../lib/x64/Release/;%LIB%
set LIB=C:\ProgramData\C\lib\;%LIB%

swig -Wall -c++ -python submodular.i
if %ErrorLevel% GTR 0 goto BadCopy

set LIB=C:\Python27\libs;C:\Python25\libs;C:\Python26\libs;%LIB%
set INCLUDE=C:\Python27\include;C:\Python25\include;C:\Python26\include;C:\Python27\Lib\site-packages\numpy\core\include;%INCLUDE%
cl /EHsc /LD /MD /Zi /O2 /D "NDEBUG" /D "WIN32" submodular.cpp submodular_wrap.cxx grd.lib thirdpartylib.lib libClp.lib libCoinUtils.lib /Fe_submodular.pyd
if %ErrorLevel% GTR 0 goto BadCopy

del *.lib
del *.obj
del *.exp
del *.pdb
del *.ilk

echo Everything done!


:BadCopy

set LIB=%OLDLIB%
set INCLUDE=%OLDINCLUDE%
