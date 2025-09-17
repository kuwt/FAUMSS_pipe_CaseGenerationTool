import sys
import os
import json
import pathlib
import shutil
import subprocess
import glob
import ast

if __name__ == "__main__":
    try:
        json_input_file_path = sys.argv[1]
        print("load json path = {} ".format(json_input_file_path))
        parametername = sys.argv[2]
        print("parametername = {} ".format(parametername))
        parametervalue = sys.argv[3]
        print("parametervalue = {} ".format(parametervalue))        
        parametertype = sys.argv[4]
        print("parametertype = {} ".format(parametertype))   
    except:
        raise Exception("no 4 arguments")

    # operation on values
    if parametertype == "bool":
        if parametervalue == "true" or parametervalue == "True" or parametervalue == "1":
            parametervalue = True
        else:
            parametervalue = False

    if parametertype == "numerics":
        parametervalue = ast.literal_eval(parametervalue)

    pathList = []
    for json_file_path in glob.glob(json_input_file_path ):
        pathList.append(json_file_path)
    for json_file_path in pathList:
        with open(json_file_path) as f:
            data = json.load(f)
                
        data[parametername] = parametervalue

        with open(json_file_path, "w") as f:
            json.dump(data, f, indent=4)
        
        print("convert {}  done".format(json_file_path))
