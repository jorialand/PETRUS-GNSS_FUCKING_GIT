#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "           CODE  RATE  DETECTION               " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

#REJECTION_CAUSE["MAX_CODE_RATE"]=9
echo "SOLUTION"
gawk '$8==9' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==9' PREPRO_OBS_TLSZ_Y19D014.dat
