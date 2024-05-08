import numpy as np
from itertools import permutations
from itertools import product
import numpy as np

characters_actual = ["A", "B", "C", "D", "E", "F"]
transition_actual = [("A", "B"), ("A", "C"), ("A", "D"), ("A", "E"), ("A", "F"), ("A", "A"), ("B", "A"),
                    ("C", "B"), ("D", "C"), ("E", "D"), ("F", "E")]

def is_valid(itinerary, transition_graph):
    for i in range(0, len(itinerary) - 1):
        if (itinerary[i], itinerary[i+1]) not in transition_graph:
            return False
    return True

def is_periodic(itinerary):
    if itinerary[0] == itinerary[len(itinerary) - 1]:
        return True
    return False

def is_repeated_substring(little_string, big_string):
    string_concat = little_string
    for i in range(1, int(len(big_string)/len(little_string)) + 1):
        if string_concat == big_string:
            return True
        else:
            string_concat += little_string
    return False


def is_shifted_itinerary(baseline, compared):
    baseline = list(baseline)
    compared = list(compared)
    for i in range(len(compared) + 1):
        if not baseline == compared:
            compared_last = str(compared[len(compared) - 1])
            for i in range(len(compared) - 1, 0, -1):
                compared[i] = compared[i - 1]
            compared[0] = compared_last
        else:
            return True
    return False

def shift_letter_to_first(letter, compared):
    compared = list(compared)
    for i in range(len(compared) + 1):
            if not compared[0] == letter:
                compared_last = str(compared[len(compared) - 1])
                for i in range(len(compared) - 1, 0, -1):
                    compared[i] = compared[i - 1]
                compared[0] = compared_last 
            else:
                new_str = ""
                for element in compared:
                    new_str += element
                return new_str
    

def generate_itineraries(period, characters, transition_graph, first_letter = None):
    if first_letter == None:
        first_letter = characters[0]
    # Keeps track of all of the valid periodic itineraries
    all_transitions = []
    

    for i in range(2, period + 2):
        # Gets all permutations of the characters as a list and converts to string
        all_possible = list(product(characters, repeat=i))
        string_itineraries = []
        for permutation in all_possible:
            itinerary_string = result = ''.join(permutation)
            string_itineraries.append(itinerary_string)
        
        # Checks if a itinerary is valid based on our transitions and periodic
        valid_transitions = []
        for it in string_itineraries:
            if is_valid(it, transition_graph) and is_periodic(it):
                valid_transitions.append(it[0:len(it)- 1])
        
        # Checks if it is directly repeated, or has a smaller period hiding inside it
        for itinerary in valid_transitions:
            if all_transitions == []:
                all_transitions.append(itinerary)
            else:
                isRepeated = False
                for smaller in all_transitions:
                    if is_repeated_substring(smaller, itinerary) or  (itinerary in all_transitions):
                        isRepeated = True
                if not isRepeated:
                    all_transitions.append(itinerary)

    # Checks if 2 itineraries are just offset versions of each other
    for baseline in all_transitions:
            for it in all_transitions:
                if (not it==baseline) and is_shifted_itinerary(it, baseline):
                    all_transitions.remove(it)
    for baseline in all_transitions:
            for it in all_transitions:
                if (not it==baseline) and is_shifted_itinerary(it, baseline):
                    all_transitions.remove(it)
                    print(f"removing({it})")

    # Shifts so that if we can, an it will have a preferred leading letter
    final_transitions = []
    for element in all_transitions:
        final_transitions.append(shift_letter_to_first(first_letter, element))
    return final_transitions
    return all_transitions


itineraries = generate_itineraries(7, characters=characters_actual, transition_graph=transition_actual)
print(itineraries)