#!/bin/bash
#
# Launch Web Interface
# Maintained by: Prabir
#

echo "╔════════════════════════════════════════════════════════╗"
echo "║                                                        ║"
echo "║     Clinical NGS Pipeline - Web Interface             ║"
echo "║              Maintained by: Prabir                     ║"
echo "║                                                        ║"
echo "╚════════════════════════════════════════════════════════╝"
echo ""

# Change to web directory
cd "$(dirname "$0")/web"

# Check Python dependencies
echo "🔍 Checking dependencies..."
python3 -c "import flask, flask_cors" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠️  Missing dependencies. Installing..."
    pip3 install flask flask-cors
fi

echo ""
echo "🚀 Starting web server..."
echo ""
echo "═══════════════════════════════════════════════════════════"
echo "  Web Interface: http://localhost:5000"
echo "  API Endpoint:  http://localhost:5000/api"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Start server
python3 api_server.py
