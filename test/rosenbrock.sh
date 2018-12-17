#!/bin/bash

if [ $# -lt 2 ]; then
    echo "usage: $0 x0 x1 iterations=10 delta=1"
    echo
    echo "easy example:  $0  1 1"
    echo "hard example:  $0 -1 1 400"
    exit -1
fi

x0=$1
x1=$2
iterations=${3:-10}
delta=${4:-1}

rm -f neldermead.state

echo "$x0 $x1" >parameters.dat

for i in $(seq $iterations)
do
    read x0 x1 <parameters.dat
    echo "(1-($x0))^2 + 100*($x1-($x0)*($x0))^2" | bc -l >energy.dat
    ../neldermead parameters.dat energy.dat $delta
done

echo
tail -v -n+0 neldermead.log neldermead.state parameters.dat min-*
