"""Create boundary conditions to test boundaryhandler

"""

import numpy as np


def writeBoundaryFile(name, string, x,y):
    f = open(name,"w")
    f.write(string)
    for i in range(len(x)):
        values = "%f  %f\n"%(x[i],y[i])
        f.write(values)
    f.close()

    return


##################### Create the boundary files ########################
info_string = "#first line 	\"bdtype\" can be one of the following: depth discharge hq-curve \
              \n#second line 	\"bdid\" gives the id of the boundary \
              \n#third line (optional?)	\"initialvalue\" \"intialh\" and \"relaxcoefficient\"  (mainly needed for state discharge curve) \
              \n#further lines  a list of two entries x y, inflow is negative, depth is given as h+z\n"

#dummy boundary file
depth_boundary_name = "outflow.dat"
depth_boundary_string = info_string + "bdtype depth \n bdid 2\n "
depth_boundary_x = [0,10.0,2400,3600]
depth_boundary_y = [0.5,0.5,1.0,1.0]
writeBoundaryFile(depth_boundary_name,depth_boundary_string,depth_boundary_x ,depth_boundary_y)

#depth_boundary_name = "outflow.dat"
#depth_boundary_string = info_string + "bdtype discharge \n bdid 2\n "
#depth_boundary_x = [0,120.0,400.0,3600]
#depth_boundary_y = [-0.5,-0.5,-0.0,0.0]
#writeBoundaryFile(depth_boundary_name,depth_boundary_string,depth_boundary_x ,depth_boundary_y)


#q boundary file
#q_boundary_name = "inflow.dat"
#q_boundary_string =  info_string + "bdtype discharge \n bdid 1\n "
#q_boundary_x = [0,120.0,400.0,3600]
#q_boundary_y = [-0.5,-0.5,-0.0,-0.0]
#writeBoundaryFile(q_boundary_name,q_boundary_string,q_boundary_x ,q_boundary_y)

q_boundary_name = "inflow.dat"
q_boundary_string =  info_string + "bdtype depth \n bdid 1\n "
q_boundary_x = [0,10.0,2400,3600]
q_boundary_y =  [0.5,0.5,1.0,1.0]
writeBoundaryFile(q_boundary_name,q_boundary_string,q_boundary_x ,q_boundary_y)


#################### Create boundary_boxes.dat #####################
boxes_string = "#definition of a boundary box with x0,x1,y0,y1,id as integer, filename \n"
#boxes_string += "199.9 201 -1.0 2000.0 2 " + depth_boundary_name +"\n"
#boxes_string += "-100.0 0.01 -1.0 2000.0 1 " +  q_boundary_name +  "\n"

f = open("boundary_boxes.dat","w")
f.write(boxes_string)
f.close()
