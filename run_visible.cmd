@echo off
setlocal
set "DISTRO=Ubuntu-24.04"
for %%A in (a b c d e f g h i j k l m n o p q r s t u v w x y z) do if /I "%~d0"=="%%A:" set "WIN_DRIVE=%%A"
set "LINUX_PROJECT=/mnt/%WIN_DRIVE%%~p0"
set "LINUX_PROJECT=%LINUX_PROJECT:\=/%"
if not defined LINUX_PROJECT (
  echo Failed to resolve the project path in WSL.
  pause
  exit /b 1
)
for /f "delims=" %%I in ('powershell.exe -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"') do set "RUN_ID=%%I"
if not defined RUN_ID (
  echo Failed to create a unique run ID.
  pause
  exit /b 1
)
echo Starting a fresh DOE in this visible terminal with twelve one-rank trials. The dashboard stays here and solver logs remain per trial.
wsl.exe -d %DISTRO% --cd "%LINUX_PROJECT%" -- python3 -u -m optimizer.calibration_controller --run --parallel-trials 12 --run-root "calibration_runs_%RUN_ID%"
set "RC=%ERRORLEVEL%"
echo.
echo Campaign process exited with code %RC%.
pause
exit /b %RC%
