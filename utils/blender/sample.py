import sys
import bpy
context = bpy.context 
scene = context.scene

cube1 = context.active_object
cube2 = cube1.copy()
# create a copy of mesh data, for bool mod
cube2.data = cube1.data.copy()
cube1.scale = (2.0, 2.0, 0.9)
mod = cube1.modifiers.new("SomeName", type='BOOLEAN')
mod.operation = 'DIFFERENCE'
# uncomment below to have cube2 in scene
#scene.objects.link(cube2)
mod.object = cube2
bpy.ops.object.modifier_apply(apply_as='DATA', modifier=mod.name)

print(sys.argv)
