import bpy
 
# Create a simple cube.
bpy.ops.mesh.primitive_cube_add()
 
# Resize the cube.
bpy.ops.transform.resize(value=(5, 3, 0.5))
 
# Get the cube object and rename it.
cube = bpy.context.object
cube.name = 'cube'
 
# Create a simple cylinder.
bpy.ops.mesh.primitive_cylinder_add(radius = 1)
 
# Get the cylinder object and rename it.
cyl = bpy.context.object
cyl.name = 'cylinder'
 
# Change the location of the cylinder.
cyl.location = (5, 3, 0)
 
# Create a boolean modifier named 'my_bool_mod' for the cube.
mod_bool = cube.modifiers.new('my_bool_mod', 'BOOLEAN')
# Set the mode of the modifier to DIFFERENCE.
mod_bool.operation = 'DIFFERENCE'
# Set the object to be used by the modifier.
mod_bool.object = cyl  
 
 
# The modifier_apply function only works on the active object.
# Set the cube as the active object.
bpy.context.scene.objects.active = cube
 
# Apply the modifier.
res = bpy.ops.object.modifier_apply(modifier = 'my_bool_mod')
 
# Delete the cylinder.
cyl.select = True
bpy.ops.object.delete()
