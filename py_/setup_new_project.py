import argparse
import os, sys


def setupDir(setup_handle, dir_bool):
    """
    make empty files and directories

    input
        setup_handle (str) - full path to directory or file to be made
        dir_bool (bool) - if True: make directory, else False: touch file
    """
    
    # check if path exists
    if os.path.exists(setup_handle) is False: 
        
        # make directory
        if dir_bool is True:
            os.mkdir(setup_handle)
        
        # make file
        elif dir_bool is False:
            os.system(f"touch {setup_handle}")

# 
def main(argv):

    # argparse
    arg_parser = argparse.ArgumentParser(description=" describe argparse")

    arg_parser.add_argument("project_name", type="str", default="MyNewProject", help='name folder for new project')

    args = arg_parser.parse_args()

    NAME = args.project_name

    # make the project name dir

    PATH = os.path.join(os.getcwd(), NAME)
    FILES = [
            ("config.ini", False)
            ("bin", True),
            ("results", True),
            ("data", True),
         ]
    for handle, dir_bool in FILES:
        setup_handle = os.path.join(PATH, handle)
        
        setupDir(setup_handle, dir_bool)
        
        
        # write to config. 


    
if __name__ == "__main__":
    main(sys.argv[1:])