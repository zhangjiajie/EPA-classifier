#! /usr/bin/env python
try:
    import sys
    import os
    from epac.ete2 import Tree, SeqGroup
    from epac.json_util import RefJsonParser, RefJsonChecker
except ImportError, e:
    print("Some packages are missing, please re-downloand EPA-classifier")
    print e
    sys.exit()

if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print "Usage: ./read_refjson JSON_FILE [field_name]"
        sys.exit()

    refjson_fname = sys.argv[1]
    if len(sys.argv) > 2:
        field_name = sys.argv[2]
    else:
        field_name = "metadata"
    
    try:
        refjson = RefJsonParser(refjson_fname)
        refjson.validate()
    except ValueError:
        print("Invalid json file format!")
        sys.exit()

    field_str = refjson.get_field_string(field_name)
    if field_str:
        print field_str
    else:
        print "Field not found: %s" % field_name
    
