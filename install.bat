@echo off
echo Installing BioX...

python --version >nul 2>&1
if errorlevel 1 (
    echo Error: Python not found. Please install Python first
    exit /b 1
)

pip --version >nul 2>&1
if errorlevel 1 (
    echo Warning: pip not found. Trying alternative installation methods...
    
    conda --version >nul 2>&1
    if not errorlevel 1 (
        echo Found conda, installing dependencies...
        conda install -y numpy>=1.19.0 tqdm>=4.45.0 multiprocess>=0.70.0
    ) else (
        echo Error: No package manager found
        echo Required packages:
        echo - numpy ^>= 1.19.0
        echo - tqdm ^>= 4.45.0
        echo - multiprocess ^>= 0.70.0
        exit /b 1
    )
) else (
    echo Installing dependencies using pip...
    pip install numpy>=1.19.0 tqdm>=4.45.0 multiprocess>=0.70.0
)

set "INSTALL_DIR=%USERPROFILE%\BioX"
if not exist "%INSTALL_DIR%" mkdir "%INSTALL_DIR%"

xcopy /E /I /Y src "%INSTALL_DIR%\src"

echo @echo off > "%INSTALL_DIR%\biox.bat"
echo python -c "import sys; sys.path.append(r'%INSTALL_DIR%'); from src.biox import main; main()" %%* >> "%INSTALL_DIR%\biox.bat"

setx PATH "%PATH%;%INSTALL_DIR%"

echo BioX installation completed!
echo Please restart your terminal, then run 'biox --help' to see usage instructions
pause 