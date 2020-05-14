#cython: wraparound=False
import numpy as np
cimport numpy as np

# Steffen, M. (1990). A simple method for monotonic interpolation in one dimension. Astronomy and Astrophysics, 239, 443.
def steffen_3d(np.ndarray[double,ndim=3, mode="c"] input_data,np.ndarray[double,ndim=3, mode="c"] input_levels,np.ndarray[double, ndim=1,mode="c"] output_level_array):
    cdef Py_ssize_t i_max=input_data.shape[0]
    cdef Py_ssize_t j_max=input_data.shape[1]
    cdef Py_ssize_t k_max=input_data.shape[2]
    cdef Py_ssize_t k_max_output=output_level_array.shape[0]
    cdef Py_ssize_t k_max_minus=k_max-1
    cdef double[:,:,:] yp=np.zeros((i_max,j_max,k_max),dtype=np.double, mode="c")
    cdef double[:,:,:] output_data=np.zeros((i_max,j_max,k_max_output),dtype=np.double, mode="c")

    cdef Py_ssize_t i,j,k,k_out,k_temp,k_high,k_low
    cdef double delta_lower,delta_upper,i_slope_lower,i_slope_upper,pp
    cdef double delta,i_slope,a,b,c,d,t,t2,t3

    for i in range(i_max):
        for j in range(j_max):
            # first point
            delta_lower = input_data[i,j,1] - input_data[i,j,0]
            delta_upper = input_data[i,j,2] - input_data[i,j,1]
            i_slope_lower = (input_levels[i,j,1] - input_levels[i,j,0])/delta_lower
            i_slope_upper = (input_levels[i,j,2] - input_levels[i,j,1])/delta_upper
            pp = i_slope_lower*(1 + delta_lower/(delta_lower + delta_upper)) - i_slope_upper*delta_lower/(delta_lower + delta_upper)
            if(pp*i_slope_lower <= 0.0):
                yp[i,j,0] = 0.0
            elif(np.abs(pp) > 2*np.abs(i_slope_lower)):
                yp[i,j,0] = 2.0*i_slope_lower
            else:
                yp[i,j,0] = pp
                
            # intermediate points
            for k in range(1,k_max_minus):
                delta_lower = input_data[i,j,k] - input_data[i,j,k-1]
                delta_upper = input_data[i,j,k+1] - input_data[i,j,k]
                i_slope_lower = (input_levels[i,j,k] - input_levels[i,j,k-1])/delta_lower        
                i_slope_upper = (input_levels[i,j,k+1] - input_levels[i,j,k])/delta_upper
                pp = (i_slope_lower*delta_upper + i_slope_upper*delta_lower)/(delta_lower + delta_upper)
               
                if (i_slope_lower*i_slope_upper <= 0.0):
                   yp[i,j,k] = 0.0
                elif (np.abs(pp) > 2.0*np.abs(i_slope_lower)):
                   yp[i,j,k] = np.copysign(2.0,i_slope_lower)*min(np.abs(i_slope_lower),np.abs(i_slope_upper))
                elif (np.abs(pp) > 2.0*np.abs(i_slope_upper)):
                   yp[i,j,k] = np.copysign(2.0,i_slope_lower)*min(np.abs(i_slope_lower),np.abs(i_slope_upper))
                else:
                   yp[i,j,k] = pp

            # last point
            delta_lower = input_data[i,j,k_max_minus-1] - input_data[i,j,k_max_minus-2]
            delta_upper = input_data[i,j,k_max_minus] - input_data[i,j,k_max_minus-1]
            i_slope_lower = (input_levels[i,j,k_max_minus-1] - input_levels[i,j,k_max_minus-2])/delta_lower
            i_slope_upper = (input_levels[i,j,k_max_minus] - input_levels[i,j,k_max_minus-1])/delta_upper
            pp = i_slope_upper*(1 + delta_upper/(delta_upper + delta_lower)) - i_slope_lower*delta_upper/(delta_upper + delta_lower)
            if (pp*i_slope_upper <= 0.0):
                yp[i,j,k_max_minus] = 0.0
            elif (np.abs(pp) > 2.0*np.abs(i_slope_upper)):
                yp[i,j,k_max_minus] = 2.0*i_slope_upper
            else:
                yp[i,j,k_max_minus] = pp

            # loop over output points
            k_temp=0
            for k_out in range(1,k_max_output):
                while((k_temp<k_max) and (input_data[i,j,k_temp]<output_level_array[k_out])):
                     k_temp=k_temp+1
                if(k_temp>0 and k_temp<k_max):
                    k_high=k_temp
                    k_low=k_high-1
                    delta = input_data[i,j,k_high] - input_data[i,j,k_low]
                    i_slope = (input_levels[i,j,k_high] - input_levels[i,j,k_low])/delta
                    a = (yp[i,j,k_low] + yp[i,j,k_high] - 2*i_slope)/(delta*delta)
                    b = (3*i_slope - 2*yp[i,j,k_low] - yp[i,j,k_high])/delta
                    c = yp[i,j,k_low]
                    d = input_levels[i,j,k_low]
                    t = output_level_array[k_out] - input_data[i,j,k_low]
                    t2 = t*t
                    t3 = t2*t
                    output_data[j,i,k_out] = a*t3 + b*t2 + c*t + d
                else:
                    output_data[j,i,k_out]= np.nan
    return output_data
