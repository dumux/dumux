import gmsh

####
#containers for the geometry entities
points = {}
lines = {}
surfaces = {}

def p(row,col):
    """return the coordinates of the point
    """
    if f"{row}{col}" not in points:
        raise ValueError("Undefined point!")
    return points[f"{row}{col}"]

def l(dir, start_row, start_col):
    """return the line object.
    dir 0 in horizontal direction,
    dir 1 in vertical direction.
    """
    key = f"{dir}{start_row}{start_col}"
    if dir > 1:
        raise ValueError("Wrong direction!")
    if key not in lines:
        raise ValueError("Undefined line!")
    return lines[key]

def s(start_row, start_col):
    """ return the surface object.
    start point is the left upper corner.
    The curve loop runs counter-clockwise.
    """
    key = f"{start_row}{start_col}"
    if key not in surfaces:
        raise ValueError("Undefined surface!")
    return surfaces[key]

def newPoint(row,col,x,y):
    """add a new point

    Args:
        row (): row index
        col (): col index
        x (): x coordinate component
        y (): y coordinate component
    """
    key = f"{row}{col}"
    if key in points:
        raise ValueError("Point already exists!")
    points[key] = gmsh.model.geo.addPoint(x,y,0)

def newLine(dir, start_row, start_col):
    """add a new line

    Args:
        row: row index of start point
        col: col index of start point
        dir: 0 for horizontal, 1 for vertical
    """
    key = f"{dir}{start_row}{start_col}"
    if key in lines:
        raise ValueError("Line exists!")
    start = p(start_row, start_col)
    end = p(start_row+1, start_col) if dir == 1 else p(start_row, start_col+1)
    lines[key] = gmsh.model.geo.addLine(start,end)

def newSurface(start_row, start_col):
    key = f"{start_row}{start_col}"
    if key in surfaces:
        raise ValueError("Surface exists!")
    l1 = l(1,start_row, start_col) # left line
    l2 = l(0,start_row+1, start_col) # bottom line
    l3 = -l(1,start_row, start_col+1) # right line
    l4 = -l(0,start_row, start_col) #top line
    curveLoop = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4])
    surfaces[key] = gmsh.model.geo.addPlaneSurface([curveLoop])

def removePoint(pointRowIdx, pointColIdx):
    pointId = p(pointRowIdx, pointColIdx)
    gmsh.model.geo.remove([(0,pointId)])

def showGeometry():
    gmsh.model.geo.synchronize()
    gmsh.fltk.run()
