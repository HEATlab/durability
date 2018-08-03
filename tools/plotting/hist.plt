#!/bin/gnuplot

# Variables to pass:
#   binwidth (float)
#   filename (string)

bin(x,width)=width*floor(x/width)

plot filename using (bin($1,binwidth)):(1.0) smooth freq with boxes


