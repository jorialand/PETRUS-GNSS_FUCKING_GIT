#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "           MASK ANGLE  REJECTION               " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["MASKANGLE"]=2
echo "SOLUTION"
gawk '$8==2' SOLUTION.dat | head

echo ""
echo "RESULT"
gawk '$8==2' PREPRO_OBS_TLSZ_Y19D014.dat | head
