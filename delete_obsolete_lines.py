#!/usr/bin/env python
# coding: utf-8

import sys

if len(sys.argv)<3:
    print("format: "+sys.argv[0]+" <input file name> <output file name>")
else:
    ifname = sys.argv[1]
    ofname = sys.argv[2]
    openstring = "# 1:T "
    closestring = "# -- restart MC: found lower energy"
    
    # find line numbers
    startlinenum = -1
    endlinenum = -1
    with open(ifname, "r") as fp:
        i=0
        for line in fp:
            if line.startswith(openstring):
                startlinenum = i
            if line.startswith(closestring):
                endlinenum = i
            i+=1
            
    #delete lines
    with open(ifname, "r") as ifp:
        with open(ofname, "w") as ofp:
            i=0
            for line in ifp:
                if i<=startlinenum or i>endlinenum:
                    ofp.write(line)
                i+=1
