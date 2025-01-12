#!/bin/bash

# Give permission: chmod +x Run.sh

# Run make to compile the code
make

# Check if the make command was successful
if [ $? -eq 0 ]; then
    # Run the executable created by make
    ./run

    # Check if the executable ran successfully
    if [ $? -eq 0 ]; then
        echo "Execution completed successfully."
    else
        echo "Execution failed."
    fi
else
    echo "Make command failed."
fi

# Clean up by running make clean
make clean
