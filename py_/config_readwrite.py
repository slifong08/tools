import configparser

def read_config(name):
    
    if ".ini" in name:
        configfile_name = name
    else:
        configfile_name = f"{name}.ini"  # name the config file

    config = configparser.ConfigParser()  # call configparser
    
    config.read(configfile_name)  # read the config file 
    
    return config, configfile_name

def write_config(config, configfile_name):  # write the config file
    
    with open(configfile_name, 'w') as configfile:  # open the configfile in write mode
        
        config.write(configfile)  # write the config object to a configfile
        
    configfile.close()  # close
    
    
def check_section(config, section):
    if config.has_section(section) is False:
        config.add_section(section)

        
def read(name):
    # short0cut for read_config
    config, configfile_name = read_config(name)
    return config, configfile_name


def write(config, configfile_name):  
    # shortcut for write the config file
    write_config(config, configfile_name)

    
def check(config, section):
    # shortcut for check_section
    check_section(config, section)

    
def checkOpt(config, section, option):
    """
    check if option is in config. 
    returns False if not in seection
    """
    exists = config.get(section, option,  fallback=False)
    
    return exists

def writeConfigDict(in_dict, config, section):
    """
    write dictionary key and values to config section
    
    input 
        in_dict (dict) - dictionary of key value pairs to write
        config (config file)
        section (str) - name of section
        
    method
        1. parse key value pairs, write to config
    
    return 
        config
    """
    
    for key, value in in_dict.items():
        config[section][key]=value

    return config

def addDict(in_dict, config, section):
    # shortcut for write configdict
    writeConfigDict(in_dict, config, section)