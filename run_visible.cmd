@echo off
setlocal
set "DISTRO=Ubuntu-24.04"
for %%A in (a b c d e f g h i j k l m n o p q r s t u v w x y z) do if /I "%~d0"=="%%A:" set "WIN_DRIVE=%%A"
set "LINUX_PROJECT=/mnt/%WIN_DRIVE%%~p0"
set "LINUX_PROJECT=%LINUX_PROJECT:\=/%"
echo NEW campaign: one sequential repose/opening case per trial; 12 trials x 1 rank.
echo Launch remains blocked until the single-cut measurement preview is approved.
wsl.exe -d %DISTRO% --cd "%LINUX_PROJECT%" -- python3 -u -m optimizer.opening_controller --run --parallel-trials 12
set "RC=%ERRORLEVEL%"
echo Campaign process exited with code %RC%.
pause
exit /b %RC%
