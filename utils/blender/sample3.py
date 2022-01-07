import bpy

context = bpy.context
scene = context.scene

bpy.ops.mesh.primitive_cube_add()
bpy.ops.mesh.primitive_cylinder_add(radius = 1)

cube = scene.objects.get("Cube")
cyl = scene.objects.get("Cylinder")
if cube and cyl:
    bool = cube.modifiers.new(name='booly', type='BOOLEAN')
    bool.object = cyl
    bool.operation = 'DIFFERENCE'
    bpy.ops.object.modifier_apply(
            {"object": cube},
            apply_as='DATA',
            modifier=bool.name)
