data_dirpath="/mnt/c/Users/mic/dev/mat"
out_dirpath="/mnt/c/Users/mic/dev/adv"

gltf2bam ${data_dirpath}/result_terrain__0_0.glb ${out_dirpath}/result_terrain__0_0.bam
bam2egg -o ${out_dirpath}/result_terrain__0_0.egg ${out_dirpath}/result_terrain__0_0.bam
egg-trans -np -t -o ${out_dirpath}/result_terrain__0_0_ed.egg ${out_dirpath}/result_terrain__0_0.egg
egg2bam -o ${out_dirpath}/result_terrain__0_0_ed.bam ${out_dirpath}/result_terrain__0_0_ed.egg

gltf2bam ${data_dirpath}/result_terrain__0_1.glb ${out_dirpath}/result_terrain__0_1.bam
bam2egg -o ${out_dirpath}/result_terrain__0_1.egg ${out_dirpath}/result_terrain__0_1.bam
egg-trans -np -o ${out_dirpath}/result_terrain__0_1_ed.egg ${out_dirpath}/result_terrain__0_1.egg
egg2bam -o ${out_dirpath}/result_terrain__0_1_ed.bam ${out_dirpath}/result_terrain__0_1_ed.egg


gltf2bam ${data_dirpath}/result_terrain__1_0.glb ${out_dirpath}/result_terrain__1_0.bam
gltf2bam ${data_dirpath}/result_terrain__1_1.glb ${out_dirpath}/result_terrain__1_1.bam

gltf2bam ${data_dirpath}/result_buildings__0_0.glb ${out_dirpath}/result_buildings__0_0.bam
gltf2bam ${data_dirpath}/result_buildings__0_1.glb ${out_dirpath}/result_buildings__0_1.bam
gltf2bam ${data_dirpath}/result_buildings__1_0.glb ${out_dirpath}/result_buildings__1_0.bam
gltf2bam ${data_dirpath}/result_buildings__1_1.glb ${out_dirpath}/result_buildings__1_1.bam