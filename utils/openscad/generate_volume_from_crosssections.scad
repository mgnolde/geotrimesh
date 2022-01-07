hollow_scale = 0.999;

dirpath="D://adv//assets//cars//pimpmydrawing";

//filename_base_suffix="bmw-x5-vector-cars-94-pimpmydrawing";
//filename_base_suffix="audi-q7-cad-cars-93-pimpmydrawing";
//filename_base_suffix="mercedes-sclass-vector-cars-112-pimpmydrawing";
//filename_base_suffix="tford-vector-cars-95-pimpmydrawing";
//filename_base_suffix="volvo-s40-vector-cars-111-pimpmydrawing";
filename_base_suffix="vw-golf-vector-cars-92-pimpmydrawing";


module extruded_side_view(filepath_side){
        color([1,0,0])
        rotate([90,0,90])
        linear_extrude(height=500)
        //import(file = "audi-q7-cad-cars-93-pimpmydrawing_side.dxf");
        import(file = filepath_side);
}


module extruded_front_view(filepath_front){

        color([0,0,1])
        rotate([90,0,0])
        linear_extrude(height=750.0)
        //import(file = "audi-q7-cad-cars-93-pimpmydrawing_front.dxf");
        import(file = filepath_front);
}    

module extruded_top_view(filepath_top){

        color([0,1,0])
        rotate([180,0,90])
        linear_extrude(height=250.0)
        //import(file = "audi-q7-cad-cars-93-pimpmydrawing_top.dxf");
        import(file = filepath_top);
}    


module solid_body(filepath_side,filepath_front,filepath_top)
{

	intersection(){	
        translate([-250,0,0]) extruded_side_view(filepath_side);
        translate([0,375,0]) extruded_front_view(filepath_front);
        translate([0,0,125]) extruded_top_view(filepath_top);
	}

}

module solid_body_test(filepath_side,filepath_front,filepath_top)
{

	union(){	
        translate([-250,0,0]) extruded_side_view(filepath_side);
        translate([0,375,0]) extruded_front_view(filepath_front);
        translate([0,0,125]) extruded_top_view(filepath_top);
	}

}


module car(filepath_side,filepath_front,filepath_top)
{

	difference(){
		solid_body(filepath_side,filepath_front,filepath_top);
		//solid_body_test(filepath_side,filepath_front,filepath_top);

        //scale([hollow_scale,hollow_scale,hollow_scale])
        //solid_body(filepath_side,filepath_front,filepath_top);
	}
}


filepath_side=str(dirpath, "//", filename_base_suffix, "_side", ".dxf");
filepath_front=str(dirpath, "//", filename_base_suffix, "_front", ".dxf");
filepath_top=str(dirpath, "//", filename_base_suffix, "_top", ".dxf");

car(filepath_side,filepath_front,filepath_top);
