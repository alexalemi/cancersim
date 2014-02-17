
#version 3.6;

#include "colors.inc"
#include "textures.inc"
#include "shapes.inc"
#include "metals.inc"
#include "glass.inc"
#include "woods.inc"


global_settings { max_trace_level 6 assumed_gamma 1.00 ambient_light 0.05}

#declare asse_z=45;
#declare sample_d=35;

camera { 
   //orthographic
    angle 35
   location < 30, 3, -110 >
   up    <0,1,0>
   right  <1,0,0>		
   look_at <30,30,0>
}


// Uncomment the area lights only if you've got lots of time.*
                   
light_source {< 30, 30, -50> color White
     fade_distance 150 fade_power 2
   area_light <40, 0, 0>, <0, 40, 0>, 30, 30
   adaptive 1
   jitter
}  


      
      
sky_sphere {
     pigment {Black}
    /*  gradient y
      color_map {
        [ 0.1  color CornflowerBlue ]
        [ 1.0  color MidnightBlue ]
      }
      scale 2
      translate -10}*/
}  
      
     


            
#declare Floor =
plane { y, -0.5
 //  pigment { color Grey*0.2 }
//   pigment {  checker LightWood, MediumWood scale 4 }   
     texture{ Brushed_Aluminum }
     finish { reflection 0.1 } 
//   texture{ Floor_Texture_1 }
}
     
     
     
//#declare filename="/home/hipe/CANCER/CancerSim/above/indaepidermis_final_state_3D.txt" 
#declare filename="input.txt"            
           
      

//#declare phongvalue=1.0;
#declare specularvalue=1.0;
//#declare phongsize=120;
#declare roughnessvalue=0.01;
#declare phongreflection=0.0; 

#declare phongvalue_cyl=0.0;
#declare phongsize_cyl=0.0;
#declare phongreflection_cyl=0.0; 
  

 
    

#declare Spheres =
union {
#include filename
      }


union {
    //object { Floor }
    object { Spheres }
}  