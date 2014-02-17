#!/bin/bash
echo Please, enter the type of initial condition 0'--->'dermis, 1'--->'epidermis, 2'--->'top of basal membrane, 3'--->'bottom of basal membrane
read init 
if [ $init == 0 ] 
then
    convert -delay 20 -loop 0 2D/indadermis_*.png  2D/indadermis_2D_animated.gif
    rm 2D/indadermis_*png
    ./3Dimages_dermis.py
    echo 'rendering '  indadermis_3D_0.txt 
    cp 3D/indadermis_3D_0.txt input.txt
    povray  making_3D_png.pov
    mv making_3D_png.png 3D/indadermis_3D_.00000.png
    for i in `seq 1 10`; do
 	echo 'rendering '  indadermis_3D_$i.txt 
  	cp 3D/indadermis_3D_$i.txt input.txt
	povray  making_3D_png.pov
 	var=$(echo "$i*0.00001" | bc)
 	mv making_3D_png.png 3D/indadermis_3D_$var.png
    done
    convert -delay 20 -loop 0 3D/indadermis_3D*.png  3D/indadermis_3D_animated.gif
    rm 3D/indadermis_3D*txt
    rm 3D/indadermis_3D*png
fi
if [ $init == 1 ] 
then
    convert -delay 20 -loop 0 2D/indaepidermis_*.png  2D/indaepidermis_2D_animated.gif
    rm 2D/indaepidermis_*png
    ./3Dimages_epidermis.py
    echo 'rendering '  indaepidermis_3D_0.txt 
    cp 3D/indaepidermis_3D_0.txt input.txt
    povray  making_3D_png.pov
    mv making_3D_png.png 3D/indaepidermis_3D_.00000.png
    for i in `seq 1 10`; do
	echo 'rendering '  indaepidermis_3D_$i.txt 
	cp 3D/indaepidermis_3D_$i.txt input.txt
	povray  making_3D_png.pov
	var=$(echo "$i*0.00001" | bc)
	mv making_3D_png.png 3D/indaepidermis_3D_$var.png
    done
    convert -delay 20 -loop 0 3D/indaepidermis_3D*.png  3D/indaepidermis_3D_animated.gif
    rm 3D/indaepidermis_3D*txt
    rm 3D/indaepidermis_3D*png
fi
if [ $init == 2 ]
then
    convert -delay 20 -loop 0 2D/top_basal_membrane_*.png  2D/top_basal_membrane_2D_animated.gif
    rm 2D/top_basal_membrane_*png
    ./3Dimages_top_basal_membrane.py
    echo 'rendering '  top_basal_membrane_3D_0.txt 
    cp 3D/top_basal_membrane_3D_0.txt input.txt
    povray  making_3D_png.pov
    mv making_3D_png.png 3D/top_basal_membrane_3D_.00000.png
    for i in `seq 1 10`; do
	echo 'rendering '  top_basal_membrane_3D_$i.txt 
	cp 3D/top_basal_membrane_3D_$i.txt input.txt
	povray  making_3D_png.pov
	var=$(echo "$i*0.00001" | bc)
	mv making_3D_png.png 3D/top_basal_membrane_3D_$var.png
    done
    convert -delay 20 -loop 0 3D/top_basal_membrane_3D*.png  3D/top_basal_membrane_animated.gif
    rm 3D/top_basal_membrane_3D*txt
    rm 3D/top_basal_membrane_3D*png
fi
if [ $init == 3 ]
then
    convert -delay 20 -loop 0 2D/bottom_basal_membrane_*.png  2D/bottom_basal_membrane_2D_animated.gif
    rm 2D/bottom_basal_membrane_*png
    ./3Dimages_bottom_basal_membrane.py
    echo 'rendering '  bottom_basal_membrane_3D_0.txt 
    cp 3D/bottom_basal_membranen_3D_0.txt input.txt
    povray  making_3D_png.pov
    mv making_3D_png.png 3D/bottom_basal_membrane_3D_.00000.png
    for i in `seq 1 10`; do
	echo 'rendering '  bottom_basal_membrane_3D_$i.txt 
	cp 3D/bottom_basal_membrane_3D_$i.txt input.txt
	povray  making_3D_png.pov
	var=$(echo "$i*0.00001" | bc)
	mv making_3D_png.png 3D/bottom_basal_membrane_3D_$var.png
    done
    convert -delay 20 -loop 0 3D/bottom_basal_membrane_3D*.png  3D/bottom_basal_membrane_animated.gif
    rm 3D/bottom_basal_membrane_3D*txt
    rm 3D/bottom_basal_membrane_3D*png
fi
rm input.txt    