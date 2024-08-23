import json

def subpeaks(path_to_json, phred_threshold = 15, subpeak_length = 5):
    """
    Extracts subpeaks from DNA basecalling data based on quality threshold and subpeak length.
    tracy basecall -f json

    Parameters:
    path_to_json (str): Path to the JSON file containing basecalling data.
    phred_threshold (int): Minimum quality threshold for considering a base as a subpeak. Defaults to 15.
    subpeak_length (int): Minimum number of consecutive subpeaks required to consider a sequence as valid. Defaults to 5.

    Returns:
    tuple: A tuple containing:
        - is_subpeaks (bool): True if there is a sequence of subpeaks that meets the length criteria, False otherwise.
        - subpeaks_positions (list): A list of tuples where each tuple contains:
            - position (str): The position in the sequence where the subpeak occurs.
            - subpeak (str): The subpeak sequence at that position.
    """
    with open(path_to_json, 'r') as file:
        data = json.load(file)

    basecallQual = data["basecallQual"]
    basecalls = data["basecalls"]

    count = 0
    is_subpeaks = False
    subpeaks_positions = []
    for item in zip(basecalls.values(), basecallQual):
        position = item[0].split(":")[0]
        base = item[0].split(":")[1]
        quality = item[1]
        #print(position, base, quality, item)
        length = len(base)
        if length > 1:
            if quality >= phred_threshold:
                count += 1
                subpeak = "|".join(base.split("|")[1:])
                subpeaks_positions.append((position, subpeak))
                #print("subpeak", subpeak, count)
                if count >= subpeak_length:
                    #print("mix", count)
                    is_subpeaks = True
                else:
                    pass
            else:
                count = 0
        else:
            count = 0
    return is_subpeaks, subpeaks_positions

        
# path = "test/results.json"
# is_subpeaks, subpeaks_positions = subpeaks(path)
# print(subpeaks_positions)




