#!/bin/bash

# Build the HTML docs into a scratch directory first. Only wipe and repopulate
# ../../docs once that build has actually succeeded, so a broken build can never
# leave the published docs empty.

make clean
mkdir -p temp_build/html

if make html; then
    rm -rf ../../docs/*
    mv -f temp_build/html/* ../../docs/
    touch ../../docs/.nojekyll
else
    echo "make html failed; leaving ../../docs untouched." >&2
    exit 1
fi
