def domain_filename_parse(date,prefix_str,directories_dict,domain_input):
    if type(domain_input) == dict:
        filename=directories_dict['output_domains_dir']+'/'+prefix_str+'_'+date+'_'+domain_input['name']+'.nc'
    if type(domain_input) == str:
        filename=directories_dict['output_domains_dir']+'/'+prefix_str+'_'+date+'_'+domain_input['name']+'.nc'
    return filename

def trajectory_filename_parse(directories_dict,trajectory_input):
    if type(trajectory_input) == dict:
        filename=directories_dict['output_trajectories_dir']+'/traj_'+trajectory_input['name']+'.nc'
    if type(trajectory_input) == str:
        filename=directories_dict['output_trajectories_dir']+'/traj_'+trajectory_input+'.nc'
    return filename
    
def forcings_filename_parse(directories_dict,forcings_input):
    if type(forcings_input) == dict:  
        filename=directories_dict['output_forcings_dir']+'/for_'+forcings_input['name']+'.nc'
    if type(trajectory_input) == str:
        filename=directories_dict['output_forcings_dir']+'/for_'+forcings_input+'.nc'
    return filename
