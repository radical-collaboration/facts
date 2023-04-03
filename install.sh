#!/bin/bash
if [[ "$VIRTUAL_ENV" != "" ]]
then
    echo "venv is active $VIRTUAL_ENV"
    echo "collecting package requirments..."
    pip freeze
else
    echo "Please activate your venv"
fi