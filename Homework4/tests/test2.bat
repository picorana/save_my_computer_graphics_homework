

@echo off

rem ..\bin\Release\pathtrace 01_textured.json

rem ..\bin\Release\pathtrace 02_area.json

rem ..\bin\Release\pathtrace 03_env.json

rem ..\bin\Release\pathtrace 04_light.json

rem ..\bin\Release\pathtrace 05_materials.json

rem ..\bin\Release\pathtrace 06_cb_direct.json

rem ..\bin\Release\pathtrace 07_cb_indirect.json

rem ..\bin\Release\pathtrace 08_test1.json


    rem ******************  MAIN CODE SECTION
    set STARTTIME=%TIME%

    rem Your code goes here (remove the ping line)
..\bin\Release\pathtrace macbook.json


    set ENDTIME=%TIME%

    rem ******************  END MAIN CODE SECTION


    rem Change formatting for the start and end times
    for /F "tokens=1-4 delims=:.," %%a in ("%STARTTIME%") do (
       set /A "start=(((%%a*60)+1%%b %% 100)*60+1%%c %% 100)*100+1%%d %% 100"
    )

    for /F "tokens=1-4 delims=:.," %%a in ("%ENDTIME%") do (
       set /A "end=(((%%a*60)+1%%b %% 100)*60+1%%c %% 100)*100+1%%d %% 100"
    )

    rem Calculate the elapsed time by subtracting values
    set /A elapsed=end-start

    rem Format the results for output
    set /A hh=elapsed/(60*60*100), rest=elapsed%%(60*60*100), mm=rest/(60*100), rest%%=60*100, ss=rest/100, cc=rest%%100
    if %hh% lss 10 set hh=0%hh%
    if %mm% lss 10 set mm=0%mm%
    if %ss% lss 10 set ss=0%ss%
    if %cc% lss 10 set cc=0%cc%

    set DURATION=%hh%:%mm%:%ss%,%cc%

    echo Start    : %STARTTIME%
    echo Finish   : %ENDTIME%
    echo          ---------------
    echo Duration : %DURATION% 