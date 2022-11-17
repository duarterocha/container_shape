#!/bin/sh
rm -rf 'bifurcation_file.txt'
cd parallel_runs
cat $(find . -name 'bifurcation_file.txt') | sort -g > bifurcation_file.txt