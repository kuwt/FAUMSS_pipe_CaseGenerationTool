import json
###########################################################
# load json parameter
############################################################
def read(json_file_path, json_key):
    try:
        with open(json_file_path) as f:
            d = json.load(f)
    except:
        raise Exception("fail to load json file {}".format(json_file_path))
    
    try:
        read_value = d[json_key]
        return read_value
    except:
        raise Exception("There is no json key ->{}<-".format(json_key))
    
def readwithdefault(json_file_path, json_key, default_value):
    try:
        with open(json_file_path) as f:
            d = json.load(f)
    except:
        raise Exception("fail to load json file {}".format(json_file_path))
    
    try:
        read_value = d[json_key]
        return read_value
    except:
        print("There is no json key ->{}<-. Use default value {}".format(json_key,default_value))
        return default_value