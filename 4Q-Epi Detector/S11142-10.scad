
difference(){

difference(){
    union(){
        translate([20,15/2-36.2/2,0.1]){
            difference(){
                translate([5.5,0,0])
                cube([11.5,36.2,1.5]);
                translate([8.5,4.4,-0.1])
                cylinder(2,d=2.2,$fn=20);
                translate([17-3,4.4,-0.1])
                cylinder(2,d=2.2,$fn=20);
                translate([8.5,36.2-4.4,-0.1])
                cylinder(2,d=2.2,$fn=20);
                translate([17-3,36.2-4.4,-0.1])
                cylinder(2,d=2.2,$fn=20);
            }
        }
    }
    translate([32.5,15/2-12/2,0])
    cube([20,12,1.7]);
}


        

union(){
    difference(){
        union(){
            color("white")
            cube([35,15,0.5]);
            translate([0.5,0.5,0.5])
            color("black")
            cube([14,14,0.25]);
            translate([7.5-0.15/2,0.5,0.75])
            color("silver")
            cube([0.15,14,0.01]);
            translate([0.5,7.5-0.15/2,0.75])
            color("silver")
            cube([14,0.15,0.01]);
            color("silver")
            translate([0,0,0.5])
            cube([16.5,15,0.01]);
        }
        translate([7.5,7.5,0.4])
        color("black")
        cylinder(2,d=2,$fn=20);
        translate([7.5,7.5,-0.1])
        color("black")
        cylinder(0.6,d=7,$fn=20);
    }
    
    translate([35-1.8-0.3,7.5-1.2/2,0.5])
    color("gold")
    cube([1.8,1.2,0.01]);
    
    translate([35-1.8-0.3,2.54+7.5-1.2/2,0.5])
    color("gold")
    cube([1.8,1.2,0.01]);
    translate([35-1.8-0.3,2*2.54+7.5-1.2/2,0.5])
    color("gold")
    cube([1.8,1.2,0.01]);
    translate([35-1.8-0.3,-2.54+7.5-1.2/2,0.5])
    color("gold")
    cube([1.8,1.2,0.01]);
    translate([35-1.8-0.3,-2*2.54+7.5-1.2/2,0.5])
    color("gold")
    cube([1.8,1.2,0.01]);
}
}