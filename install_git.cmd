@echo off
where git >nul 2>&1
if %ERRORLEVEL%==0 (
    exit /b 0
)

echo Git not found. Detecting package manager...

where winget >nul 2>&1
if %ERRORLEVEL%==0 (
    echo Installing Git...
    winget install --id Git.Git -e --accept-package-agreements --accept-source-agreements
    if %ERRORLEVEL% neq 0 (
        echo winget failed to install Git.
        exit /b 3
    )
    goto :verify
)

where choco >nul 2>&1
if %ERRORLEVEL%==0 (
    echo Installing Git...
    choco install git -y
    if %ERRORLEVEL% neq 0 (
        echo choco failed to install Git.
        exit /b 4
    )
    goto :verify
)

where scoop >nul 2>&1
if %ERRORLEVEL%==0 (
    echo Installing Git...
    scoop install git
    if %ERRORLEVEL% neq 0 (
        echo scoop failed to install Git.
        exit /b 5
    )
    goto :verify
)

echo No supported package manager found (winget, choco, scoop).
echo Please install Git manually: https://git-scm.com/download/win
exit /b 1

:verify
where git >nul 2>&1
if %ERRORLEVEL%==0 (
    echo Installation complete.
    exit /b 0
) else (
    echo Installation completed but git was not found in PATH. Check installer output.
    exit /b 2
)