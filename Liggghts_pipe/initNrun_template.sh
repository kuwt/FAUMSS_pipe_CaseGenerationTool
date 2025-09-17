#!/bin/bash

echo "checking if init is previously done"

if [ -f "./restart/restart.bk1" ]; then
    echo "initialization is done before. Skip initialization."
else
    echo "init"
    lmp_auto -in init.in
fi

if [ -f "./restart/restart.bkfinal" ]; then
	echo "mainrun is done before. Skip mainrun."
else
	echo "main run"
	mpirun -np {numOfProcessor} lmp_auto -in run.in
fi

echo "detail run"
mpirun -np {numOfProcessor} lmp_auto -in detailrun.in
