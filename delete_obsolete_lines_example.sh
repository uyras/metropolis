#!/bin/bash
# this example removes obsolete lines from 
# all files in folder "before"
# and saves result in folder "after"
for f in before/*.out; do ./removeUnusedResults.py $f
"after/"$(basename $f); done