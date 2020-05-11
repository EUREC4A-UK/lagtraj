def get_paths():
    with open('../../default_directories') as input_file:
        paths_dict=yaml.load(this_case_file, Loader=yaml.FullLoader)
    
