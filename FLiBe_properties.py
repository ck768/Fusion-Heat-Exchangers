import numpy as np
import festim as F
import h_transport_materials as htm

import numpy as np
import matplotlib.pyplot as plt


def calculate_inlet_velocity(
    flow_rate,
    inlet_diameter,
    breeder_density,
    breeder,
    volumetric=False,
    suppress_print=True,
):
    """Calculate the inlet velocity of fluid breeder at a given flow rate, inlet diameter, breeder density, and temperature.
    Used for OpenFOAM simulation.

    Parameters
    ----------
    flow_rate : float
        Flow rate in kg/s.
    inlet_diameter : float
        Inlet diameter in m.
    breder_density : float
        Breeder fluid density in kg/m3.
    breeder : str
        Breeder fluid name.
    volumetric : bool
        True if flow rate is volumetric. False if flow rate is mass flow rate.
    suppress_print : bool
        Supresses print output when True.

    Returns
    -------
    float
        Inlet velocity in m/s.
    """

    inlet_area = np.pi * (inlet_diameter / 2) ** 2  # m^2

    if volumetric:
        inlet_velocity = flow_rate / inlet_area
    else:
        inlet_velocity = flow_rate * breeder_density ** (-1) * inlet_area ** (-1)  # m/s
    if not suppress_print:
        print(
            f"Inlet velocity for {breeder} flow rate of {flow_rate}kg/s is {inlet_velocity}m/s."
        )

    return inlet_velocity


def calculate_reynolds_number(
    inlet_velocity,
    characteristic_length,
    kinematic_viscosity,
    breeder,
    suppress_print=True,
):
    """Calculate the reynolds number of a fluid breeder given its inlet velocity, characteristic length, and kinematic viscosity.
    Used to determine if turbulence is needed for OpenFOAM simulation.

    Parameters
    ----------
    inlet_velocity : float
        Inlet velocity in m/s.
    characteristic_length : float
        Characteristic length in m.
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
     breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Reynolds number (dimensionless).
    """

    reynolds_number = (inlet_velocity * characteristic_length) / kinematic_viscosity

    if not suppress_print:
        print(f"Reynolds number for {breeder} is {reynolds_number}.")

        if reynolds_number > 3500:
            print(f"Flow is turbulent.")
        else:
            print(f"Flow is laminar.")

    return reynolds_number


def calculate_schmidt_number(
    kinematic_viscosity,
    diffusivity,
    breeder,
    suppress_print=True,
):
    """Calculate the schmidt number of a fluid breeder given its kinematic viscosity and mass diffusivity.
    Used for turbulent diffusion term in FESTIM.

    Parameters
    ----------
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
    diffusivity : float
        Fluid diffusivity in m2/s.
    breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Schmidt number (dimensionless).
    """

    schmidt_number = kinematic_viscosity / diffusivity

    if not suppress_print:
        print(f"Schmidt number for {breeder} is {schmidt_number}.")

    return schmidt_number


def plot_reynolds_number_vs_inlet_velocity(
    characteristic_length,
    kinematic_viscosity,
    breeder_temperature,
    breeder,
    model_velocity,
    show=True,
):
    """Plot the reynolds number of a fluid breeder with varying inlet velocities.
    Assumes constant characteristic length and kinematic viscosity.

    Parameters
    ----------
    characteristic_length : float
        Characteristic length in m.
    kinematic_viscosity : float
        Kinematic viscosity in m2/s.
    breeder_temperature : float
        Breeder temperature in K.
     breeder : str
        Breeder fluid name.
    """
    inlet_velocities = np.linspace(0, 2, 100000)  # m/s
    Re_numbers = []

    for inlet_velocity in inlet_velocities:  # m/s
        Re_numbers.append(
            calculate_reynolds_number(
                inlet_velocity,
                characteristic_length,
                kinematic_viscosity,
                breeder,
                suppress_print=True,
            )
        )

    plt.plot(inlet_velocities, Re_numbers, "b-")
    plt.axhline(y=3500, color="r", linestyle="--", label="Turbulence Threshold")
    plt.axvline(x=model_velocity, color="g", label="Model Inlet Velocity")
    plt.xlabel("Inlet Velocity (m/s)")
    plt.ylabel("Reynolds Number")
    plt.title(
        f"Reynolds Number vs Inlet Velocity for {breeder} at {breeder_temperature}K."
    )
    plt.legend()
    if show:
        plt.show()


# reference: https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-epsilon.html


def calculate_initial_k(inlet_velocity, suppress_print=True):
    """Calculate initial turbulence kinetic energy.

    Parameters
    ----------
    inlet_velocity : float
        Inlet velocity in m/s.

    Returns
    -------
    float
        Initial turbulence kinetic energy in m2/s2.
    """

    # assume initial turbulence is isotropic so that U'2_x = U'2_y = U'2_z
    k = (3 / 2) * (0.05 * inlet_velocity) ** 2  # 5% of inlet velocity

    if not suppress_print:
        print(f"Initial turbulence kinetic energy for FLiBe: {k} m2/s2")

    return k


def calculate_initial_epsilon(k, characteristic_length, suppress_print=True):
    """Calculate initial turbulence dissipation rate.

    Parameters
    ----------
    k : float
        Initial turbulence kinetic energy in m2/s2.
    characteristic_length : float
        Characteristic length in m.

    Returns
    -------
    float
        Initial turbulence dissipation rate in m2/s3.
    """

    epsilon = 0.09 ** (3 / 4) * k ** (3 / 2) / characteristic_length

    if not suppress_print:
        print(f"Initial turbulence dissipation rate for FLiBe: {epsilon} m2/s3")

    return epsilon


def calculate_initial_omega(k, characteristic_length, suppress_print=True):
    """Calculate initial specific dissipation rate.

    Parameters
    ----------
    k : float
        Initial turbulence kinetic energy in m2/s2.
    characteristic_length : float
        Characteristic length in m.

    Returns
    -------
    float
        Initial specific dissipation rate in 1/s.
    """

    omega = np.sqrt(k) / (0.09 ** (1 / 4) * characteristic_length)

    if not suppress_print:
        print(f"Initial specific dissipation rate for FLiBe: {omega} m2/s3")
    return omega


def calculate_FLiBe_kinematic_viscosity(
    breeder_temperature,
    breeder_density,
    breeder,
    suppress_print=True,
):
    """Calculate the kinematic viscosity of FLiBe at a given temperature and density.
    Used for OpenFOAM simulation.

    Parameters
    ----------
    breeder_temperature : float
        Breeder temperature in K.
    breeder_density : float
        Breeder density in kg/m3.
     breeder : str
        Breeder fluid name.

    Returns
    -------
    float
        Kinematic viscosity in m2/s.
    """
    breeder_dynamic_viscosity = 0.000116 * np.exp(
        3755 / breeder_temperature
    )  # Pa.s = kg s / m s2; equation from Williams 2006

    kinematic_viscosity = breeder_dynamic_viscosity / breeder_density  # m2/s

    if not suppress_print:
        print(
            f"Kinematic viscosity of {breeder} at {breeder_temperature}K is {kinematic_viscosity}m2/s."
        )

    return kinematic_viscosity


breeder = "FLiBe"

breeder_temperature = 908  # K from Meschini 2021
FLiBe_density = 2245 - 0.424 * (
    breeder_temperature - 273.15
)  # kg/m3 ; equation from Vidrio 2022

print(FLiBe_density)

flow_rate = 1.5 / 3600 # m3/h -> to m3/s from https://www.sciencedirect.com/science/article/pii/S1359431116303738

inlet_diameter = 0.03 # m from CAD

k_b = F.k_B  # eV/K, boltzmann constant

flibe_diffusivity = htm.diffusivities.filter(material=htm.FLIBE).mean()
E_D = flibe_diffusivity.act_energy.magnitude  # eV
D_0 = flibe_diffusivity.pre_exp.magnitude  # m2/s

FLiBe_diffusivity = D_0 * np.exp(-E_D / (k_b * breeder_temperature))  # m2/s

inlet_velocity = calculate_inlet_velocity(
    flow_rate, inlet_diameter, FLiBe_density, breeder, volumetric=True
)
print(inlet_velocity)

kinematic_viscosity = calculate_FLiBe_kinematic_viscosity(
    breeder_temperature, FLiBe_density, breeder
)

Re = calculate_reynolds_number(
    inlet_velocity, inlet_diameter, kinematic_viscosity, breeder, suppress_print=False
)

k = calculate_initial_k(inlet_velocity, suppress_print=False)
epsilon = calculate_initial_epsilon(k, characteristic_length=inlet_diameter)
omega = calculate_initial_omega(
    k, characteristic_length=inlet_diameter, suppress_print=False
)

plot_reynolds_number_vs_inlet_velocity(
    inlet_diameter, kinematic_viscosity, breeder_temperature, breeder, inlet_velocity, show=False
)

dynamic_viscosity = 6e-3
specific_heat = 2400
thermal_conductivity = 1

show_print = False
prandtl_number = dynamic_viscosity * specific_heat / thermal_conductivity
if show_print == True:
    print(f"prandtl number is {prandtl_number}")
