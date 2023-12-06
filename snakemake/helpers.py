def return_parsed_extra_params(extra_params):
    #starting blank
    cmd = ""
    #for key in extra parameters
    for key in extra_params:
        #append the key value pair if it's a parmeter that needed something
        if extra_params[key]:
            cmd += " --{0} {1}".format(key,extra_params[key])
        else: #otherwise if it's parameter that's just a flag, append just the flag
            cmd += " --{0}".format(key)
    return(cmd)
