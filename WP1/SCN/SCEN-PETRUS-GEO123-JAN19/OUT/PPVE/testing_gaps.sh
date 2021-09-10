#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "              GAPS DETECTION                   " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["DATA_GAP"]=6
echo "SOLUTION"
gawk '$8==6' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==6' PREPRO_OBS_TLSZ_Y19D014.dat
