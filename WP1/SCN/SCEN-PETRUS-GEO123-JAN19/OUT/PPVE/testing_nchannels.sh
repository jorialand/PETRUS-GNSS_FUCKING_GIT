#!/bin/bash

echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "             TESTING WP1-PPVE                  " 
echo "         MAX NCHANNELS REJECTION               " 
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo ""

head -1 SOLUTION.dat 

# REJECTION_CAUSE["NCHANNELS_GPS"]=1
echo "SOLUTION"
gawk '$8==1' SOLUTION.dat

echo ""
echo "RESULT"
gawk '$8==1' PREPRO_OBS_TLSZ_Y19D014.dat
