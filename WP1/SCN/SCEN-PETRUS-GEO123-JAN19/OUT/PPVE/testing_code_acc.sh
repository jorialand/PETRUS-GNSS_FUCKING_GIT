#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "            CODE ACC  DETECTION                " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["MAX_CODE_RATE_STEP"]=10
echo "SOLUTION"
gawk '$8==10' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==10' PREPRO_OBS_TLSZ_Y19D014.dat
