
    // Parameters
    thickness = 10; // Plate thickness
    holeDiameter = 5; // Hole diameter

    // Module to create the main plate with holes
    module mainPlate() {
        difference() {
            cube([100, 100, thickness], true); // Create the base plate
            // Drill holes
            for (pos = [[10, 10], [20, 20]]) {
                translate([pos[0], pos[1], -thickness/2])
                    cylinder(h = thickness * 2, d = holeDiameter, $fn = 100); // Make holes through the plate
            }
        }
    }

    // Call the main module to create the plate
    mainPlate();
    