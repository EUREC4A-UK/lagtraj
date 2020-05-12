def domain_filename_parse(date,prefix_str,directories_dict,domain_dict):
    filename=directories_dict['output_domains_dir']+'/'+prefix_str+'_'+date+'_'+domain_dict['name']+'.nc'
    return filename

def trajectory_filename_parse(directories_dict,trajectory_dict):
    filename=directories_dict['output_trajectories_dir']+'/traj_'+trajectory_dict['name']+'.nc'
    return filename
    
def forcing_filename_parse(directories_dict,forcing_dict):
    filename=directories_dict['output_forcings_dir']+'/for_'+forcing_dict['name']+'.nc'
    return filename
