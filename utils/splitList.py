def split_list(input_list, num_splits):
    avg = len(input_list) / float(num_splits)
    output = {}
    last = 0.0
    index = 0
    while last < len(input_list):
        output[index] = input_list[int(last):int(last + avg)]
        last += avg
        index += 1
    return output