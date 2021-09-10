#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "            MAX PSR REJECTION                  " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["MAX_PSR_OUTRNG"]=4
echo "SOLUTION"
gawk '$8==4' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==4' PREPRO_OBS_TLSZ_Y19D014.dat
