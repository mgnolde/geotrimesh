#!/bin/bash

#DATA_DIRPATH="/mnt/e/devenv/geotrimesh/data"
DATA_DIRPATH="/media/hd/devenv/geotrimesh/data"
PROJECT_NAME="zh"
SRS="EPSG:2056"


TILES_VERTICAL="0 1 2 3 4 5 6 7 8 9 "
TILES_HORIZONTAL="0 1 2 3 4 5 6 7 8 9"
TILE_HEIGHT=160
TILE_WIDTH=160


mkdir -p ${DATA_DIRPATH}/proc/${PROJECT_NAME}/feat
mkdir -p ${DATA_DIRPATH}/proc/${PROJECT_NAME}/adm
mkdir -p ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho
mkdir -p ${DATA_DIRPATH}/proc/${PROJECT_NAME}/dem
mkdir -p ${DATA_DIRPATH}/res/${PROJECT_NAME}


STEPS=""
#STEPS=${STEPS}"acquire_data",
#STEPS=${STEPS}"get_dependencies",
STEPS=${STEPS}"perform_preprocessing",
STEPS=${STEPS}"perform_processing"
#STEPS=${STEPS}"perform_postprocessing"
#STEPS="all"


STEP="get_dependencies"
if [[ ${STEPS} == *${STEP}* ]] || [[ ${STEPS} = "all" ]]
then

    echo "[STEP] "${STEP}
    sudo apt install bc jq imagemagick assimp-utils meshlab

fi


STEP="acquire_data"
if [[ ${STEPS} == *${STEP}* ]] || [[ ${STEPS} = "all" ]]
then

    echo "[STEP] "${STEP}

    #LoD2 Building data
    #https://www.stadt-zuerich.ch/geodaten/download/3D_Dachmodell_LoD2

fi


STEP="perform_preprocessing"
if [[ ${STEPS} == *${STEP}* ]] || [[ ${STEPS} = "all" ]]
then


    ## convert buildings data from dxf to stl. Since ASSIMP seems to export invalid STL files,
    ## Meshlab as used for that step
     
    #assimp \
    #  export \
    #  ${DATA_DIRPATH}/orig/${PROJECT_NAME}/feat/buildings/data.dxf \
    #  ${DATA_DIRPATH}/proc/${PROJECT_NAME}/feat/buildings/buildings.obj

    #meshlabserver \
    #  -i ${DATA_DIRPATH}/proc/${PROJECT_NAME}/feat/buildings/buildings.obj \
    #  -o ${DATA_DIRPATH}/proc/${PROJECT_NAME}/feat/buildings/buildings.stl \
    #  -m vt vn wt \


    for SOURCE in "dem" "ortho"
    do

        ## Merge and clip orthofotos

        rm ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic.vrt
        gdalbuildvrt \
          ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic.vrt \
          ${DATA_DIRPATH}/orig/${PROJECT_NAME}/${SOURCE}/*.tif
    
        rm ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic.tif
        gdal_translate \
          -stats \
          -of GTiff \
          -a_nodata -3.40282e+38 \
          -co "COMPRESS=LZW" \
          -co BIGTIFF=YES \
          ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic.vrt \
          ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic.tif
    
    
        rm ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic_clip.tif
        gdalwarp \
          -cutline ${DATA_DIRPATH}/orig/${PROJECT_NAME}/adm/dtm_26830_12460_bbox.gpkg \
          -crop_to_cutline \
          -dstalpha \
          -s_srs ${SRS} \
          ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic.tif \
          ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic_clip.tif
    
    done


fi


STEP="perform_processing"
if [[ ${STEPS} == *${STEP}* ]] || [[ ${STEPS} = "all" ]]
then

    ## must be a divider of TILE_HEIGHT and TILE_WIDTH
    TARGET_DEM_RESOLUTION=0.5
    TARGET_ORTHO_RESOLUTION=0.1
    #SOURCE_ZMIN=`gdalinfo -json -mm ${DATA_DIRPATH}/orig/${PROJECT_NAME}/dem/dtm_26830_12460.tif | jq .bands[0].computedMin`
    #SOURCE_ZMAX=`gdalinfo -json -mm ${DATA_DIRPATH}/orig/${PROJECT_NAME}/dem/dtm_26830_12460.tif | jq .bands[0].computedMax`
    #TARGET_XMIN=2682000.0
    #TARGET_YMIN=1246000.0
    #TARGET_XMAX=2684000.0
    #TARGET_YMAX=1248000.0

    SOURCE_ZMIN=404.55
    SOURCE_ZMAX=440.36
    TARGET_XMIN=2683000.0
    TARGET_XMAX=1246000.0
    TARGET_YMIN=2684000.0
    TARGET_YMAX=1247000.0

    echo ${SOURCE_ZMIN} ${SOURCE_ZMAX}



    for TILE_VERTICAL in ${TILES_VERTICAL}
    do

        for TILE_HORIZONTAL in ${TILES_HORIZONTAL}
        do

            echo ${TILE_VERTICAL} / ${TILE_HORIZONTAL} 


            TILE_HORIZONTAL_PLUS1=$((TILE_HORIZONTAL+1))
            TILE_VERTICAL_PLUS1=$((TILE_VERTICAL+1))
    
            TILE_XMIN=$(echo "scale=10; $TARGET_XMIN+$(($TILE_HORIZONTAL*$TILE_WIDTH))" | bc) 
            TILE_DEM_XMAX=$(echo "scale=10; $TARGET_XMIN+$TARGET_DEM_RESOLUTION+$(($TILE_HORIZONTAL_PLUS1*$TILE_WIDTH))" | bc) 
            TILE_ORTHO_XMAX=$(echo "scale=10; $TARGET_XMIN+$(($TILE_HORIZONTAL_PLUS1*$TILE_WIDTH))" | bc) 
    
            TILE_DEM_YMIN=$(echo "scale=10; $TARGET_YMAX-$TARGET_DEM_RESOLUTION-$(($TILE_VERTICAL_PLUS1*$TILE_HEIGHT))" | bc) 
            TILE_ORTHO_YMIN=$(echo "scale=10; $TARGET_YMAX-$(($TILE_VERTICAL_PLUS1*$TILE_HEIGHT))" | bc) 
            TILE_YMAX=$(echo "scale=10; $TARGET_YMAX-$(($TILE_VERTICAL*$TILE_HEIGHT))" | bc) 
    
    
            #TILE_XMIN_BUF=$(echo "scale=10; $TILE_XMIN-0.1" | bc) 
            #TILE_YMAX_BUF=$(echo "scale=10; $TILE_YMAX+0.1" | bc) 
            #TILE_ORTHO_XMAX_BUF=$(echo "scale=10; $TILE_ORTHO_XMAX+0.1" | bc) 
            #TILE_ORTHO_YMIN_BUF=$(echo "scale=10; $TILE_ORTHO_YMIN-0.1" | bc) 
    
    
            echo ${TILE_XMIN} ${TILE_DEM_XMAX}
            echo ${TILE_DEM_YMIN} ${TILE_YMAX}
        

    
            for SOURCE in "dem" "ortho"
            do

                if [[ ${SOURCE} = "dem" ]]
                then
                    TILE_XMAX=${TILE_DEM_XMAX}
                    TILE_YMIN=${TILE_DEM_YMIN}
                    TARGET_RESOLUTION=${TARGET_DEM_RESOLUTION}
                elif [[ ${SOURCE} = "ortho" ]]
                then
                    TILE_XMAX=${TILE_ORTHO_XMAX}
                    TILE_YMIN=${TILE_ORTHO_YMIN}
                    TARGET_RESOLUTION=${TARGET_ORTHO_RESOLUTION}
                fi

                rm ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.tif
                gdalwarp \
                  -of GTiff \
                  -co "COMPRESS=LZW" \
                  -te ${TILE_XMIN} ${TILE_YMIN} ${TILE_XMAX} ${TILE_YMAX} \
                  -tr ${TARGET_RESOLUTION} ${TARGET_RESOLUTION} \
                  -co BIGTIFF=YES \
                  -t_srs ${SRS} \
                  ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic_clip.tif \
                  ${DATA_DIRPATH}/proc/${PROJECT_NAME}/${SOURCE}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.tif

            done


            rm ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}_mirrored.png
            gdal_translate \
              -b 1 \
              -b 2 \
              -b 3 \
              -of PNG \
              -ot Byte \
              ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.tif \
              ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}_mirrored.png


            rm ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.png
            convert \
              -flip \
              ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}_mirrored.png \
              ${DATA_DIRPATH}/proc/${PROJECT_NAME}/ortho/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.png


            rm ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.scad
            python3 geotrimesh/generate_terrain.py \
              --dem ${DATA_DIRPATH}/proc/${PROJECT_NAME}/dem/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.tif \
              --zmin ${SOURCE_ZMIN} \
              --zmax ${SOURCE_ZMAX} \
              --clippoly ${DATA_DIRPATH}/orig/${PROJECT_NAME}/adm/dtm_26830_12460_bbox.gpkg \
              --out ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.scad \

            #openscad \
            #  ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.scad \
            #  -o ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.stl



            rm ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.scad
            python3 geotrimesh/generate_terrain.py \
              --dem ${DATA_DIRPATH}/proc/${PROJECT_NAME}/dem/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}.tif \
              --zmin ${SOURCE_ZMIN} \
              --zmax ${SOURCE_ZMAX} \
              --dem_extrude_down 500 \
              --dem_extrude_up 500 \
              --clippoly ${DATA_DIRPATH}/orig/${PROJECT_NAME}/adm/dtm_26830_12460_bbox.gpkg \
              --out ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}_bbox.scad \
              -b

            #openscad \
            #  ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}_bbox.scad \
            #  -o ${DATA_DIRPATH}/res/${PROJECT_NAME}/mosaic_clip_${TILE_VERTICAL}_${TILE_HORIZONTAL}_bbox.stl



        done

    done

fi

