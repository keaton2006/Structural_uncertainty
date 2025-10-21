#!/bin/sh
for i in {150..180}; do
    echo $i
	python train.py --time_index $i &
done