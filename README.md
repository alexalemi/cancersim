# nevi_code_plos1

Code for an elastic simulation of skin cancer.

The main code is located in the `code` folder.

There are 4 executable files corresponding to each of the 4 possible initial configurations, i.e. the  position of the first nevus cell:

 * epidermis - `Sim_epidermis.py`
 * dermis - `Sim_dermis.py`
 * top of basal membrane - `Sim_top_basal_membrane.py`
 * bottom of basal membrane - `Sim_bottom_basal_membrane.py`

The code will generate as output:

 * the `.png` files inside the dir `2D/`: a graphical 2D view of the  system configuration once the stable condition has been achieved;
 * the file pressure.dat inside the dir `pressure/`: the pressure acting on each cell;
 * the file evolution.dat inside the dir `evolution/`: the position of each new cell added to the nevus;
 * the files state inside the dir `states/`: the system configuration once the stable condition has been achieved.

Once the code has finished, make sure you have downloaded the 3.6 version of [`povray`](http://www.povray.org/). Then run the script: `script_making_movies.sh`

this will generate as output 2 `gif` files:

 * a `2D.gif` inside the dir `2D/`: an animation of the nevus evolution; it will delete all the preexisting .png inside the 2D dir
 * a `3D.gif` inside the dir `3D/`: an animation of the nevus evolution rendered in 3D

