#!/bin/bash

echo "checking if init is previously done"

if [ -f "./restart/restart.bk1" ]; then
    echo "initialization is done before. Skip initialization."
else
    echo "init"
    lmp_auto -in init.in
fi

echo "main run"
mpirun -np {numOfProcessor} lmp_auto -in run.in
