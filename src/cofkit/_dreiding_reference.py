from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DreidingAtomParameters:
    r1: float
    theta0: float
    r0: float
    d0: float
    phi: float
    s: float


# Pinned from lammps-interface commit 255f027cb76142d39c050a6810404debc6a06562.
# Source: lammps_interface/dreiding.py
DREIDING_REFERENCE_SOURCE = (
    "lammps-interface commit 255f027cb76142d39c050a6810404debc6a06562 "
    "(https://github.com/peteboyd/lammps_interface)"
)


DREIDING_PARAMETERS: dict[str, DreidingAtomParameters] = {
    "H_": DreidingAtomParameters(0.33, 180.0, 3.195, 0.0152, 0.0, 12.382),
    "H__HB": DreidingAtomParameters(0.33, 180.0, 3.195, 0.0001, 0.0, 12.0),
    "H__b": DreidingAtomParameters(0.510, 90.0, 3.195, 0.0152, 0.0, 12.382),
    "B_3": DreidingAtomParameters(0.880, 109.471, 4.02, 0.095, 0.0, 14.23),
    "B_2": DreidingAtomParameters(0.790, 120.0, 4.02, 0.095, 0.0, 14.23),
    "C_3": DreidingAtomParameters(0.770, 109.471, 3.8983, 0.0951, 0.0, 14.034),
    "C_R": DreidingAtomParameters(0.700, 120.0, 3.8983, 0.0951, 0.0, 14.034),
    "C_2": DreidingAtomParameters(0.670, 120.0, 3.8983, 0.0951, 0.0, 14.034),
    "C_1": DreidingAtomParameters(0.602, 180.0, 3.8983, 0.0951, 0.0, 14.034),
    "N_3": DreidingAtomParameters(0.702, 106.7, 3.6621, 0.0774, 0.0, 13.843),
    "N_R": DreidingAtomParameters(0.650, 120.0, 3.6621, 0.0774, 0.0, 13.843),
    "N_2": DreidingAtomParameters(0.615, 120.0, 3.6621, 0.0744, 0.0, 13.843),
    "N_1": DreidingAtomParameters(0.556, 180.0, 3.6621, 0.0744, 0.0, 13.843),
    "O_3": DreidingAtomParameters(0.660, 104.51, 3.4046, 0.0957, 0.0, 13.483),
    "O_R": DreidingAtomParameters(0.660, 120.0, 3.4046, 0.0957, 0.0, 13.483),
    "O_2": DreidingAtomParameters(0.560, 120.0, 3.4046, 0.0957, 0.0, 13.483),
    "O_1": DreidingAtomParameters(0.528, 180.0, 3.4046, 0.0957, 0.0, 13.483),
    "F_": DreidingAtomParameters(0.611, 180.0, 3.4720, 0.0725, 0.0, 14.444),
    "Al3": DreidingAtomParameters(1.047, 109.471, 4.39, 0.31, 0.0, 12.0),
    "Si3": DreidingAtomParameters(0.937, 109.471, 4.27, 0.31, 0.0, 12.0),
    "P_3": DreidingAtomParameters(0.890, 93.3, 4.15, 0.32, 0.0, 12.0),
    "S_3": DreidingAtomParameters(1.040, 92.1, 4.03, 0.344, 0.0, 12.0),
    "Cl": DreidingAtomParameters(0.997, 180.0, 3.9503, 0.2833, 0.0, 13.861),
    "Ga3": DreidingAtomParameters(1.210, 109.471, 4.39, 0.4, 0.0, 12.0),
    "Ge3": DreidingAtomParameters(1.210, 109.471, 4.27, 0.4, 0.0, 12.0),
    "As3": DreidingAtomParameters(1.210, 92.1, 4.15, 0.41, 0.0, 12.0),
    "Se3": DreidingAtomParameters(1.210, 90.6, 4.03, 0.43, 0.0, 12.0),
    "Br": DreidingAtomParameters(1.167, 180.0, 3.95, 0.37, 0.0, 12.0),
    "In3": DreidingAtomParameters(1.390, 109.471, 4.59, 0.55, 0.0, 12.0),
    "Sn3": DreidingAtomParameters(1.373, 109.471, 4.47, 0.55, 0.0, 12.0),
    "Sb3": DreidingAtomParameters(1.432, 91.6, 4.35, 0.55, 0.0, 12.0),
    "Te3": DreidingAtomParameters(1.280, 90.3, 4.23, 0.57, 0.0, 12.0),
    "I_": DreidingAtomParameters(1.360, 180.0, 4.15, 0.51, 0.0, 12.0),
    "Na": DreidingAtomParameters(1.860, 90.0, 3.144, 0.5, 0.0, 12.0),
    "Ca": DreidingAtomParameters(1.940, 90.0, 3.472, 0.05, 0.0, 12.0),
    "Fe": DreidingAtomParameters(1.285, 90.0, 4.54, 0.055, 0.0, 12.0),
    "Zn": DreidingAtomParameters(1.330, 109.471, 4.54, 0.055, 0.0, 12.0),
    # These three entries are explicitly heuristic in the upstream source.
    "Cu": DreidingAtomParameters(1.302, 90.0, 4.54, 0.055, 0.0, 12.0),
    "Ni": DreidingAtomParameters(1.164, 90.0, 4.54, 0.055, 0.0, 12.0),
    "Mg": DreidingAtomParameters(1.421, 90.0, 4.54, 0.055, 0.0, 12.0),
    "C_R1": DreidingAtomParameters(0.700, 120.0, 4.23, 0.1356, 0.0, 14.034),
    "C_34": DreidingAtomParameters(0.770, 109.471, 4.2370, 0.3016, 0.0, 12.0),
    "C_33": DreidingAtomParameters(0.770, 109.471, 4.1524, 0.25, 0.0, 12.0),
    "C_32": DreidingAtomParameters(0.770, 109.471, 4.0677, 0.1984, 0.0, 12.0),
    "C_31": DreidingAtomParameters(0.770, 109.471, 3.983, 0.1467, 54.74, 12.0),
}


DREIDING_FRAMEWORK_TYPE_BY_ELEMENT: dict[str, str] = {
    "H": "H_",
    "B": "B_2",
    "C": "C_R",
    "N": "N_R",
    "O": "O_3",
    "F": "F_",
    "P": "P_3",
    "S": "S_3",
    "Cl": "Cl",
    "Br": "Br",
}

