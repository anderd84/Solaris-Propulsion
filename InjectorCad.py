import subprocess
def injector_cad_write():
    # Python code to generate an OpenSCAD script for an injector plate
    # Example parameters for the injector plate
    thickness = 10  # Plate thickness
    hole_diameter = 5  # Diameter of the holes
    hole_positions = [(10, 10), (20, 20)]  # Positions of the holes

    # Path to the OpenSCAD script file
    file_path = ('injector_plate.scad')

    # Generate the OpenSCAD script
    scad_script = f"""
    // Parameters
    thickness = {thickness}; // Plate thickness
    holeDiameter = {hole_diameter}; // Hole diameter

    // Module to create the main plate with holes
    module mainPlate() {{
        difference() {{
            cube([100, 100, thickness], true); // Create the base plate
            // Drill holes
            for (pos = [{', '.join(f'[{x}, {y}]' for x, y in hole_positions)}]) {{
                translate([pos[0], pos[1], -thickness/2])
                    cylinder(h = thickness * 2, d = holeDiameter, $fn = 100); // Make holes through the plate
            }}
        }}
    }}

    // Call the main module to create the plate
    mainPlate();
    """

    # Write the script to a file
    with open(file_path, 'w') as file:
        file.write(scad_script)

    print(f"OpenSCAD script written to {file_path}")
    
    open_command = f'openscad {file_path}'

    # Execute the command
    subprocess.call(open_command, shell=True)

    # OpenSCAD command - adjust the path to the OpenSCAD executable if needed
    open_command = f'"C:\\Program Files\\OpenSCAD\\openscad.exe" {file_path}'  # Use the full path to openscad.exe
    
    # Execute the command to open the .scad file with OpenSCAD
    subprocess.call(open_command, shell=True)