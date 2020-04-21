from typing import Tuple
import textdistance

def get_distance_multi(distance_getter: textdistance.algorithms, position_table, sequence_id_1, sequence_id_2) -> Tuple:
    ret_dict = {}

    if sequence_id_1 == sequence_id_2:
        result = 0
    else:
        result = distance_getter.distance(list(position_table.loc[sequence_id_1]),
                                                                         list(position_table.loc[sequence_id_2]))

    return (sequence_id_1, sequence_id_2, result)