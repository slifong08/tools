import configparser

def read_config(name):
    
    configfile_name = f"{name}.ini"  # name the config file

    config = configparser.ConfigParser()  # call configparser
    
    config.read(configfile_name)  # read the config file 
    
    return config, configfile_name

def write_config(config, configfile_name):  # write the config file
    
    with open(configfile_name, 'w') as configfile:  # open the configfile in write mode
        
        config.write(configfile)  # write the config object to a configfile
        
    configfile.close()  # close