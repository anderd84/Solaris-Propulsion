# =============================================================================
# Created on Nov 18, 2024
# Developing hoop stress calculations for the chamber with cooling channels
# =============================================================================

# Custom imports
from General.units import Q_, unitReg

# General imports
import numpy as np
from dataclasses import dataclass
from icecream import ic

@dataclass
class HoopStress:

    def __init__(self, chamber_pressure, cooling_pressure, outer_radius, wall_thickness, channel_diameter, num_channels):
        self.chamber_pressure = chamber_pressure
        self.cooling_pressure = cooling_pressure
        self.outer_radius = outer_radius
        self.wall_thickness = wall_thickness
        self.channel_diameter = channel_diameter
        self.num_channels = num_channels

    def effective_thickness(self):
        """
        Calculates the effective thickness after accounting for cooling channels.
        """
        self.effective_wall_thickness = self.wall_thickness - (self.channel_diameter * self.num_channels / (2 * np.pi * self.outer_radius))
        if self.effective_wall_thickness <= 0:
            raise ValueError("Effective wall thickness is non-physical (<= 0). Check your inputs.")
        return self.effective_wall_thickness


    def solid_hoop_stress(self):
        """
        Calculates the hoop stress for a solid cylinder (no cooling channels).
        """
        inner_radius = self.outer_radius - self.wall_thickness
        Dm = (self.outer_radius + inner_radius) / 2
        self.hoop_stress_solid = self.chamber_pressure * Dm / self.wall_thickness
        return self.hoop_stress_solid

    def channel_hoop_stress(self):
        """
        Calculates the hoop stress including the effects of cooling channels.
        """
        t_eff = self.effective_thickness()
        inner_radius = self.outer_radius - t_eff
        Dm = (self.outer_radius + inner_radius) / 2
        stress_from_chamber = self.chamber_pressure * Dm / t_eff
        stress_from_cooling = self.cooling_pressure * self.outer_radius / t_eff
        self.hoop_stress_channels = stress_from_chamber + stress_from_cooling
        return self.hoop_stress_channels

if __name__ == '__main__':
    # Input values with units
    chamber_pressure = Q_(500, unitReg.psi)  # Chamber pressure in psi
    cooling_pressure = Q_(100, unitReg.psi)  # Cooling pressure in psi
    outer_radius = Q_(7.75, unitReg.inch)    # Outer radius in inches
    wall_thickness = Q_(0.125, unitReg.inch) # Wall thickness in inches
    channel_diameter = Q_(0.05, unitReg.inch) # Cooling channel diameter in inches
    num_channels = 10                        # Number of cooling channels

    # Instantiate HoopStress object
    hoop_stress = HoopStress(chamber_pressure, cooling_pressure, outer_radius, wall_thickness, channel_diameter, num_channels)

    # Perform calculations
    effective_thickness = hoop_stress.effective_thickness()
    solid_stress = hoop_stress.solid_hoop_stress()
    channel_stress = hoop_stress.channel_hoop_stress()

    # Output results
    ic(effective_thickness)
    ic(solid_stress)
    ic(channel_stress)
