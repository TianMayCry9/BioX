#!/bin/bash

echo "Installing BioX..."

if ! command -v python3 &> /dev/null; then
    echo "Error: Python3 not found. Please install Python3 first"
    exit 1
fi

if ! command -v pip3 &> /dev/null; then
    echo "Warning: pip3 not found. Trying alternative installation methods..."
    
    if command -v conda &> /dev/null; then
        echo "Found conda, installing dependencies..."
        conda install -y numpy>=1.19.0 tqdm>=4.45.0 multiprocess>=0.70.0
    else
        echo "Conda not found, trying system package manager..."
        
        if [[ "$OSTYPE" == "darwin"* ]]; then
            brew install python3-numpy python3-tqdm
        elif [[ -f /etc/debian_version ]]; then
            sudo apt-get update
            sudo apt-get install -y python3-numpy python3-tqdm
        elif [[ -f /etc/redhat-release ]]; then
            sudo yum install -y python3-numpy python3-tqdm
        else
            echo "Error: Unsupported operating system"
            exit 1
        fi
    fi
else
    echo "Installing dependencies using pip..."
    pip3 install numpy>=1.19.0 tqdm>=4.45.0 multiprocess>=0.70.0
fi

INSTALL_DIR="/usr/local/bin"
if [[ ! -w "$INSTALL_DIR" ]]; then
    INSTALL_DIR="$HOME/.local/bin"
    mkdir -p "$INSTALL_DIR"
fi

BIOX_DIR="$INSTALL_DIR/biox"
mkdir -p "$BIOX_DIR"
cp -r src/* "$BIOX_DIR/"

cat > "$INSTALL_DIR/biox" << 'EOF'
#!/bin/bash
python3 -c "import sys; sys.path.append('$BIOX_DIR'); from biox import main; main()"
EOF

chmod +x "$INSTALL_DIR/biox"

echo "BioX installation completed!"
echo "Run 'biox --help' to see usage instructions" 