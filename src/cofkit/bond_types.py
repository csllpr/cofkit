from __future__ import annotations


AROMATIC_BOND_ORDER = 1.5


def normalize_bond_order(order: float) -> float:
    if abs(order - AROMATIC_BOND_ORDER) <= 0.15:
        return AROMATIC_BOND_ORDER
    rounded = round(order)
    if abs(order - rounded) <= 0.15:
        return float(rounded)
    return float(order)


def is_aromatic_bond_order(order: float) -> bool:
    return abs(normalize_bond_order(order) - AROMATIC_BOND_ORDER) <= 1.0e-6


def bond_order_to_cif_type(order: float) -> str:
    normalized = normalize_bond_order(order)
    if is_aromatic_bond_order(normalized):
        return "A"
    if abs(normalized - 1.0) <= 1.0e-6:
        return "S"
    if abs(normalized - 2.0) <= 1.0e-6:
        return "D"
    if abs(normalized - 3.0) <= 1.0e-6:
        return "T"
    return f"{normalized:.3f}".rstrip("0").rstrip(".")


def cif_type_to_bond_order(value: str | None) -> float | None:
    if value is None:
        return None
    token = value.strip()
    if not token or token in {".", "?"}:
        return None
    canonical = token.upper()
    mapping = {
        "S": 1.0,
        "SING": 1.0,
        "SINGLE": 1.0,
        "1": 1.0,
        "D": 2.0,
        "DB": 2.0,
        "DOUBLE": 2.0,
        "2": 2.0,
        "T": 3.0,
        "TRIPLE": 3.0,
        "3": 3.0,
        "A": AROMATIC_BOND_ORDER,
        "AR": AROMATIC_BOND_ORDER,
        "AROM": AROMATIC_BOND_ORDER,
        "AROMATIC": AROMATIC_BOND_ORDER,
        "DEL": AROMATIC_BOND_ORDER,
        "DELOC": AROMATIC_BOND_ORDER,
    }
    if canonical in mapping:
        return mapping[canonical]
    try:
        return normalize_bond_order(float(token))
    except ValueError:
        return None
