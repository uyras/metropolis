#!/bin/bash
eval "$@"
while [[ $? = 100 ]]; do 
  eval "$@ -f gs.mfsys"
done