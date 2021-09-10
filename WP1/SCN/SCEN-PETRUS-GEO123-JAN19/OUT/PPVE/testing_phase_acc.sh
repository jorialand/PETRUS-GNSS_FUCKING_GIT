#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "        PHASE RATE ACC  DETECTION              " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["MAX_PHASE_RATE_STEP"]=8
echo "SOLUTION"
gawk '$8==8' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==8' PREPRO_OBS_TLSZ_Y19D014.dat
