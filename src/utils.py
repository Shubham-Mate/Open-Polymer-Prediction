import re
from typing import Tuple


def get_first_number(input_string: str) -> Tuple[int, int]:
    extracted_number = -1
    continuation_index = -1

    regular_expression = r"%(\d+)"

    matches = re.search(regular_expression, input_string)

    first_match = matches.group(1)
    extracted_number = int(first_match[1:])
    continuation_index = len(first_match) - 1

    return extracted_number, continuation_index
