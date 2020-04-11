#!/bin/bash

cc='_contour'
outmap='_labelmap'
border='_border'

while read line
do
   I="${line%.*}"
   ./SphSPS -i $2$line -k $3 -m $4 -outm $I$outmap.png -outb $I$border.png -c $2$I$cc.png \n
done < $1
