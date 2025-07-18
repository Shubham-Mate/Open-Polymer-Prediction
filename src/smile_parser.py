import json
from dataclasses import dataclass
from typing import Tuple, List, Dict
from config import load_config
from smiles_constants import *


def load_chemical_data(data_path) -> Dict:
    try:
        with open(data_path, "r") as data_file:
            chemical_data = json.load(data_file)
            return chemical_data
    except:
        print("Unable to find the JSON file for chemical data")


chemical_data = load_chemical_data("./data.json")


@dataclass
class Node:
    element: str
    id: int = 0
    is_part_of_aromatic_chain: bool = False
    charge: int = 0

    def __repr__(self):
        return f"Node(id={self.id}, element='{self.element}') "


@dataclass
class Edge:
    element_id_1: int
    element_id_2: int
    bond_type: str


class Molecule:
    def __init__(self, smiles: str, data: Dict = chemical_data):
        self.smiles = smiles
        self.data = data
        self.elements = set(data.keys())
        self.nodes, self.edges = self.parse_smiles(smiles)

    def parse_smiles(self, smiles: str) -> Tuple[List[Node], List[List[Edge]]]:
        nodes = []
        edges = []
        element_count = 0
        for i in range(len(smiles)):
            char = smiles[i]
            if char.upper() in self.elements:
                curr_possible_bond_count = self.data[char.upper()]["possible_num_bonds"]
                nodes.append(Node(id=element_count, element=char.upper()))
                element_count += 1
                if i < len(smiles) - 1:
                    if smiles[i + 1] in BOND_TYPES:
                        if smiles[i + 1] == Bonds.DOUBLE_BOND.value:
                            curr_possible_bond_count -= 2
                        elif smiles[i + 1] == Bonds.TRIPLE_BOND.value:
                            curr_possible_bond_count -= 3
                    else:
                        curr_possible_bond_count -= 1

                if i > 0:
                    if smiles[i - 1] in BOND_TYPES:
                        if smiles[i - 1] == Bonds.DOUBLE_BOND.value:
                            curr_possible_bond_count -= 2
                        elif smiles[i - 1] == Bonds.TRIPLE_BOND.value:
                            curr_possible_bond_count -= 3
                    else:
                        curr_possible_bond_count -= 1
                for j in range(curr_possible_bond_count):
                    nodes.append(Node(element="H", id=element_count))
                    element_count += 1

        return nodes, edges


new_mol = Molecule("N")
print(new_mol.nodes)
