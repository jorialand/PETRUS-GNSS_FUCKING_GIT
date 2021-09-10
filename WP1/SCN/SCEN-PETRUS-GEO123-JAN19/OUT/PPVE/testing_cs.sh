#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "          CYCLE SLIPS  DETECTION               " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["CYCLE_SLIP"]=5
echo "SOLUTION"
gawk '$8==5' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==5' PREPRO_OBS_TLSZ_Y19D014.dat
