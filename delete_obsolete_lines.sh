#!/bin/bash
perl -p0e -i 's/(# 1:T[^\\n]+\\n).+# -- restart MC: found[^\\n]+\\n/$1/s' $1