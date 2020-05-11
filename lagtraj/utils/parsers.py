def domain_filename_parse(date,prefix_str,directories_dict,domain_dict):
    filename=directories_dict['output_domains_dir']+'/ml_'+date+'_'+str(domain_dict['lat_max'])+'_'+str(domain_dict['lon_min'])+'_'+str(domain_dict['lat_min'])+'_'+str(domain_dict['lon_max'])+'.nc'
    return filename
