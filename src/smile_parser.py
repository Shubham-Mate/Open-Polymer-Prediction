import json
from dataclasses import dataclass
from typing import Tuple, List, Dict
from config import load_config
from smiles_constants import *
from utils import get_first_number


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
        self.edge_adj_matrix = self.get_edge_adj_matrix()

    def parse_smiles(self, smiles: str) -> Tuple[List[Node], List[Edge]]:
        nodes = []
        edges = []

        element_ids_list = [
            None for i in range(len(smiles))
        ]  # For tracking element id for each char in smiles. If it doesn't have an id, it'll be None
        cycle_ids_mapping = {}  # A dictionary which has key as the number used to indicate connection that form cycle, and value is the list of ids of first and last element

        element_count = 0
        i = 0
        while i < len(smiles):
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

                element_ids_list[i] = element_count
                current_char_id = element_count
                element_count += 1
                items_to_be_checked = []

                # Check elements after current element to judge how many bonds of current element are used up
                if i < len(smiles) - 1:
                    items_to_be_checked.append(i + 1)
                    while smiles[
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

                # Check elements before current element to judge how many bonds of current element are used up
                if i > 0:
                    items_to_be_checked.append(i - 1)
                    if smiles[items_to_be_checked[-1]] == ")":
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

                # Final calculation of how many remaining bonds the current element can make, along with adding edges representing bonds used by other elements
                for item in items_to_be_checked:
                    if smiles[item] in BOND_TYPES:
                        curr_possible_bond_count -= BOND_TYPES.index(smiles[item]) + 1
                        bond_type = smiles[item]
                        element_id_of_prev_element = element_ids_list[item - 1]
                        while element_id_of_prev_element == None:
                            item -= 1
                            element_id_of_prev_element = element_ids_list[item - 1]
                    else:
                        curr_possible_bond_count -= 1
                        element_id_of_prev_element = element_ids_list[item]
                        bond_type = Bonds.SINGLE_BOND.value

                    if item < i:
                        edges.append(
                            Edge(
                                element_id_1=current_char_id,
                                element_id_2=element_id_of_prev_element,
                                bond_type=bond_type,
                            )
                        )

                # The remaining number of bonds current element can make, will all be with hydrogens
                for j in range(curr_possible_bond_count):
                    nodes.append(
                        Node(
                            element="H",
                            id=element_count,
                            is_part_of_aromatic_chain=char.islower(),
                        )
                    )
                    edges.append(
                        Edge(
                            element_id_1=current_char_id,
                            element_id_2=element_count,
                            bond_type=Bonds.SINGLE_BOND.value,
                        )
                    )
                    element_count += 1

            # If we come across a beginning or ending of cycle
            elif char.isnumeric() or char == "%":
                extracted_number = char

                if char == "%":
                    extracted_number, continuation_index = get_first_number(smiles[i:])
                    extracted_number = str(extracted_number)
                    i += continuation_index

                if extracted_number not in cycle_ids_mapping:
                    cycle_ids_mapping[extracted_number] = []

                j = i
                while element_ids_list[j] == None:
                    j -= 1
                cycle_ids_mapping[extracted_number].append(element_ids_list[j])

            i += 1

        for pairs in cycle_ids_mapping.values():
            edges.append(
                Edge(
                    element_id_1=pairs[0],
                    element_id_2=pairs[1],
                    bond_type=Bonds.SINGLE_BOND.value,
                )
            )

        return nodes, edges

    def get_edge_adj_matrix(self) -> List[List[Edge | None]]:
        edge_adj_matrix = [
            [None for _ in range(len(self.nodes))] for _ in range(len(self.nodes))
        ]

        for edge in self.edges:
            id_1 = edge.element_id_1
            id_2 = edge.element_id_2

            try:
                edge_adj_matrix[id_1][id_2] = edge
                edge_adj_matrix[id_2][id_1] = edge
            except:
                print(edge)

        return edge_adj_matrix


new_mol = Molecule("*CC(*)c1ccccc1C(=O)OCCCCCC")
print(new_mol.nodes)
print(new_mol.edges)
print("\n\n\n")
print(new_mol.edge_adj_matrix)
