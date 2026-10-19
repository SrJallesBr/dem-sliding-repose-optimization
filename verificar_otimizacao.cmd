@echo off
setlocal
set "DISTRO=Ubuntu-24.04"
for %%A in (a b c d e f g h i j k l m n o p q r s t u v w x y z) do if /I "%~d0"=="%%A:" set "WIN_DRIVE=%%A"
set "LINUX_PROJECT=/mnt/%WIN_DRIVE%%~p0"
set "LINUX_PROJECT=%LINUX_PROJECT:\=/%"
echo PREFLIGHT ONLY - no solver and no campaign.
wsl.exe -d %DISTRO% --cd "%LINUX_PROJECT%" -- python3 -B -m optimizer.opening_controller --preflight --parallel-trials 12
set "RC=%ERRORLEVEL%"
pause
exit /b %RC%
