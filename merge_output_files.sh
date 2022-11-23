#!/bin/sh
rm -rf 'bifurcation_file.txt'
cd parallel_runs_cylinder
cat $(find . -name 'bifurcation_file.txt') | sort -g > bifurcation_file.txt