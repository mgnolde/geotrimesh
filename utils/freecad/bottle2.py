import sys
import os

import Part, math
from FreeCAD import Base
import Mesh
import MeshPart

sys.path.append(os.path.join(os.sep, 'usr', 'lib', 'freecad', 'lib'))


App.newDocument("Unnamed")


App.setActiveDocument("Unnamed")
App.ActiveDocument=App.getDocument("Unnamed")


def makeBottleTut(myWidth = 50.0, myHeight = 70.0, myThickness = 30.0):
    aPnt1=Base.Vector(-myWidth / 2., 0, 0)
    aPnt2=Base.Vector(-myWidth / 2., -myThickness / 4., 0)
    aPnt3=Base.Vector(0, -myThickness / 2., 0)
    aPnt4=Base.Vector(myWidth / 2., -myThickness / 4., 0)
    aPnt5=Base.Vector(myWidth / 2., 0, 0)
    aArcOfCircle = Part.Arc(aPnt2, aPnt3, aPnt4)
    aSegment1=Part.LineSegment(aPnt1, aPnt2)
    aSegment2=Part.LineSegment(aPnt4, aPnt5)
    aEdge1=aSegment1.toShape()
    aEdge2=aArcOfCircle.toShape()
    aEdge3=aSegment2.toShape()
    aWire=Part.Wire([aEdge1, aEdge2, aEdge3])
    aTrsf=Base.Matrix()
    aTrsf.rotateZ(math.pi) # rotate around the z-axis
    aMirroredWire=aWire.copy()
    aMirroredWire.transformShape(aTrsf)
    myWireProfile=Part.Wire([aWire, aMirroredWire])
    myFaceProfile=Part.Face(myWireProfile)
    aPrismVec=Base.Vector(0, 0, myHeight)
    myBody=myFaceProfile.extrude(aPrismVec)
    myBody=myBody.makeFillet(myThickness / 12.0, myBody.Edges)
    neckLocation=Base.Vector(0, 0, myHeight)
    neckNormal=Base.Vector(0, 0, 1)
    myNeckRadius = myThickness / 4.
    myNeckHeight = myHeight / 10.
    myNeck = Part.makeCylinder(myNeckRadius, myNeckHeight, neckLocation, neckNormal)
    myBody = myBody.fuse(myNeck)
    return myBody
 
 


el = makeBottleTut()
Part.show(el)
el.exportStep(os.path.join(os.getcwd(), 'test.stp'))


"""
App.ActiveDocument.ActiveObject.Label = "Cylinder"
App.ActiveDocument.recompute()
 
__objs__=[]
__objs__.append(FreeCAD.getDocument("Unnamed").getObject("Cylinder"))
Mesh.export(__objs__,u"opillar_2.stl")
"""

#o = FreeCAD.getDocument("Unnamed").findObjects()[0]

#msh = MeshPart.meshFromShape(Shape=el, MaxLength=1)
#s = Part.Shape()
#s.makeShapeFromMesh(msh)
#s.exportStep()

#import FreeCAD
#import Part
#import Mesh
#shape = Part.Shape()
#shape.read('test.stp')

doc = App.newDocument('Doc')
pf = doc.addObject("Part::Feature","MyShape")
pf.Shape = el
Mesh.export([pf], 'test.stl')

#el.exportStl(os.path.join(os.getcwd(), 'test.stl'))

#del __objs__
