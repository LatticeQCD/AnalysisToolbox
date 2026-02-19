#!/bin/bash


./clean.sh

mkdir -p temp_build/html
make html

mv temp_build/html/* ../../docs/ -f
touch ../../docs/.nojekyll
