
## https://forum.freecadweb.org/viewtopic.php?style=3&t=38218

import FreeCAD
import Part
FreeCAD.newDocument("test")
FreeCAD.setActiveDocument("test")
FreeCAD.ActiveDocument = FreeCAD.getDocument("test")
# NURBS surface from a Part.BSplineSurface().
# len(knot_u) := nNodes_u + degree_u + 1
# len(knot_v) := nNodes_v + degree_v + 1
#poles=[[[15.5, 4.8, 1.46, 1.0], [10.18, 4.8, 0.89, 1.0], [1.8, 8.9, 0.57, 1.0], [0, 15, 0.43, 1.0], [0, 22, -0.02, 1.0]], [[15.15, 17.05, 1.13, 1.0], [11.7, 17.05, 0.81, 1.0], [6.72, 17.1, 0.48, 1.0], [4.31, 18.74, 0.38, 1.0], [3.23, 22.8, 0.05, 1.0]]]
poles=[[[15.5, 4.8, 1.46, 10.0], [10.18, 4.8, 0.89, 10.0], [1.8, 8.9, 0.57, 2.0], [0, 15, 0.43, 2.0], [0, 22, -0.02, 1.0]], [[15.15, 17.05, 1.13, 1.0], [11.7, 17.05, 0.81, 1.0], [6.72, 17.1, 0.48, 1.0], [4.31, 18.74, 0.38, 1.0], [3.23, 22.8, 0.05, 1.0]]]
degree_u=1
degree_v=3
nNodes_u=2
nNodes_v=5
knot_v=[0, 0.5256099762388414, 1]
knot_v_mults=[4, 1, 4]
knot_u=[0,1]
knot_u_mults=[2, 2]
NURBS_Cubic_surf=Part.BSplineSurface()
NURBS_Cubic_surf.increaseDegree(degree_u,degree_v)

for i in range(0,len(knot_u)):    #-1):
	NURBS_Cubic_surf.insertUKnot(knot_u[i],knot_u_mults[i],0.0000001)

for i in range(0,len(knot_v)):    #-1):
	NURBS_Cubic_surf.insertVKnot(knot_v[i],knot_v_mults[i],0.0000001)

poles_xyz=[]
for ii in range(0, len(poles)):
	poles_xyz.append([])
	for jj in range(0, len(poles[ii])):
		poles_xyz[ii].append(FreeCAD.Base.Vector(poles[ii][jj][0], poles[ii][jj][1], poles[ii][jj][2]))
		NURBS_Cubic_surf.setPole(ii+1, jj+1, poles_xyz[ii][jj], poles[ii][jj][3])

Part.show(NURBS_Cubic_surf.toShape())
