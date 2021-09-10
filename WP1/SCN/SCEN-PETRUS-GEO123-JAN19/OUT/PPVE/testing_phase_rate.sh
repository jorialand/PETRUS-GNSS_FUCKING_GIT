#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "           PHASE RATE  DETECTION               " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["MAX_PHASE_RATE"]=7
echo "SOLUTION"
gawk '$8==7' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==7' PREPRO_OBS_TLSZ_Y19D014.dat
