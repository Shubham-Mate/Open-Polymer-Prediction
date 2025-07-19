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
                nodes.append(
                    Node(
                        id=element_count,
                        element=char.upper(),
                        is_part_of_aromatic_chain=char.islower(),
                    )
                )
                element_count += 1
                items_to_be_checked = []
                if i < len(smiles) - 1:
                    items_to_be_checked.append(i + 1)
                    if smiles[
                        items_to_be_checked[-1]
                    ].isnumeric() and items_to_be_checked[-1] + 1 < len(smiles):
                        items_to_be_checked.append(items_to_be_checked[-1] + 1)
                    if smiles[items_to_be_checked[-1]] == "(":
                        items_to_be_checked[-1] += 1
                        open_bracket_count = 1
                        j = items_to_be_checked[-1]
                        while open_bracket_count > 0 and j < len(smiles):
                            if smiles[j] == ")":
                                open_bracket_count -= 1
                            elif smiles[j] == "(":
                                open_bracket_count += 1
                            j += 1
                        if j < len(smiles):
                            items_to_be_checked.append(j)
                    elif smiles[items_to_be_checked[-1]] == ")":
                        items_to_be_checked.pop()

                if i > 0:
                    items_to_be_checked.append(i - 1)
                    if smiles[i - 1] == ")":
                        items_to_be_checked.pop()
                        close_bracket_count = 1
                        j = i - 2
                        while close_bracket_count > 0 and j > 0 and smiles[j] != ")":
                            if smiles[j] == "(":
                                close_bracket_count -= 1
                            elif smiles[j] == ")":
                                close_bracket_count += 1
                            j -= 1
                        if j >= 0:
                            items_to_be_checked.append(j)
                    elif smiles[items_to_be_checked[-1]] == "(":
                        while smiles[items_to_be_checked[-1]] == "(":
                            items_to_be_checked[-1] -= 1

                for item in items_to_be_checked:
                    print(item, items_to_be_checked)
                    if smiles[item] in BOND_TYPES:
                        curr_possible_bond_count -= BOND_TYPES.index(smiles[item]) + 1
                    else:
                        curr_possible_bond_count -= 1
                print(curr_possible_bond_count)
                for j in range(curr_possible_bond_count):
                    nodes.append(
                        Node(
                            element="H",
                            id=element_count,
                            is_part_of_aromatic_chain=char.islower(),
                        )
                    )
                    element_count += 1

        return nodes, edges


new_mol = Molecule("CC1(CCCCC1)O")
print(new_mol.nodes)
