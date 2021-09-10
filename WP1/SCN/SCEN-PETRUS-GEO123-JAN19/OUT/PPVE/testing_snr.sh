#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "            MIN SNR REJECTION                  " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["MIN_CNR"]=3
echo "SOLUTION"
gawk '$8==3' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==3' PREPRO_OBS_TLSZ_Y19D014.dat
