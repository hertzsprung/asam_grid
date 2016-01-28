#!/bin/bash
if [[ $# -eq 0 ]] ; then
    echo 'Usage: ./archive.sh <version>'
    exit 1
fi
git archive --prefix asam-gridgen-$1/ -o ../asam-gridgen_$1.orig.tar.gz $1
