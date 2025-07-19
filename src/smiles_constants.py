from enum import Enum

BOND_TYPES = ["-", "=", "#"]


class Bonds(Enum):
    SINGLE_BOND = "-"
    DOUBLE_BOND = "="
    TRIPLE_BOND = "#"
