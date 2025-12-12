#!/bin/bash
# save as monitor_memory.sh

if [ -z "$1" ]; then
    echo "Usage: ./monitor_memory.sh <program_name>"
    echo "Example: ./monitor_memory.sh implicit_learning"
    exit 1
fi

PROGRAM_NAME=$1

echo "Waiting for $PROGRAM_NAME to start..."
while ! pgrep -x "$PROGRAM_NAME" > /dev/null; do
    sleep 1
done

PID=$(pgrep -x "$PROGRAM_NAME")
echo "Found PID: $PID"
echo "Monitoring memory usage (press Ctrl+C to stop)..."
echo ""

while kill -0 $PID 2>/dev/null; do
    MEM_KB=$(ps -p $PID -o rss= 2>/dev/null)
    if [ -n "$MEM_KB" ]; then
        MEM_GB=$(echo "scale=2; $MEM_KB/1024/1024" | bc)
        TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
        echo "[$TIMESTAMP] Memory: $MEM_GB GB ($MEM_KB KB)"
    fi
    sleep 2  # Check every 2 seconds
done

echo ""
echo "Process $PROGRAM_NAME (PID: $PID) has terminated."
