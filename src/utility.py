import numpy as np
from scipy.integrate import quad, odeint
from data import *
from figure import *

def transform(x_list,y_list,x_coordinate):
    """
    Transform the discretized PDF to continuous PDF, which makes it easier to integrate
    
    :param x_list: crystal size (um)
    :type x_list: list
    :param y_list: density function (1/um^4)
    :type y_list: list
    :param x_coordinate: x coordinate (um)
    :type x_coordinate: float
    :return: y coordinate (1/um^4)
    :rtype: float
    """
    
    y_coordinate=0.0
    if x_coordinate < x_list[0]:
        y_coordinate=y_list[0]
    elif x_coordinate > x_list[-1]:
        y_coordinate=y_list[-1]
    for i in range(len(x_list)-1):
        if x_coordinate == x_list[i]:
            y_coordinate=y_list[i]
        elif x_coordinate > x_list[i] and x_coordinate < x_list[i+1]:
            y_coordinate=y_list[i]*(x_list[i+1]-x_coordinate)/(x_list[i+1]-x_list[i])+y_list[i+1]*(x_coordinate-x_list[i])/(x_list[i+1]-x_list[i])
        elif x_coordinate == x_list[i+1]:
            y_coordinate=y_list[i+1]
    
    return y_coordinate

def moment_calculation(size,density_function,order):
    """
    Calculate the desired order moment of the specified CSD
    
    :param size: crystal size (um)
    :type size: list
    :param density_function: density function (1/um^4)
    :type density_function: list
    :param order: moment order
    :type order: float
    :return: moment
    :rtype: float 
    """
    
    x_list=size
    y_list=density_function
    continuous_function=lambda x: transform(x_list,y_list,x)*x**order
    output=quad(continuous_function,size[0],size[-1])[0]
    
    return output

def DF_to_NF(size,density_function,mesh_size):
    """
    Transform the number density function to number fraction
    
    :param size: size
    :type size: variable
    :param density_function: number density function
    :type density_function: variable
    :param mesh_size: mesh size of the size distribution (um)
    :type mesh_size: float
    :return: size and number fraction
    :rtype: variable
    """
    
    size_output=variable("size","$\mu m$",[])
    number_fraction=variable("number fraction","-",[])
    
    mesh_number=int((size.value[-1]-size.value[0])/mesh_size)+1
    tmp1=np.linspace(size.value[0],size.value[-1],mesh_number)
    normalization_constant=moment_calculation(size.value,density_function.value,0)
    
    for i in range(len(tmp1)-1):
        L_mean=(tmp1[i]+tmp1[i+1])/2
        f1=transform(size.value,density_function.value,tmp1[i])
        f2=transform(size.value,density_function.value,tmp1[i+1])
        tmp2=(f1+f2)/2*mesh_size/normalization_constant
        size_output.value.append(L_mean)
        number_fraction.value.append(tmp2)
        
    return size_output, number_fraction

def NF_to_VF(size,number_fraction):
    """
    Transform the number fraction to volume fraction
    
    :param size: size
    :type size: variable
    :param number_fraction: number fraction
    :type number_fraction: variable
    :return: size and volume fraction
    :rtype: variable
    """
    
    volume_fraction=variable("volume fraction","-",[])
    
    total_volume=0
    for i in range(len(size.value)):
        tmp1=size.value[i]**3*number_fraction.value[i]
        total_volume=total_volume+tmp1
    for i in range(len(size.value)):
        tmp1=size.value[i]**3*number_fraction.value[i]
        volume_fraction.value.append(tmp1/total_volume)    
    
    return size,volume_fraction

def VF_to_NF(size,volume_fraction):
    """
    Transform the number fraction to volume fraction
    
    :param size: size
    :type size: variable
    :param volume_fraction: volume fraction
    :type volume_fraction: variable
    :return: size and number fraction
    :rtype: variable
    """
    
    number_fraction=variable("number fraction","-",[])
    
    total_number=0
    for i in range(len(size.value)):
        tmp1=volume_fraction.value[i]/size.value[i]**3
        total_number=total_number+tmp1
    for i in range(len(size.value)):
        tmp1=volume_fraction.value[i]/size.value[i]**3
        number_fraction.value.append(tmp1/total_number)    
    
    return size,number_fraction

def array_to_list(array1):
    """
    Transform an one-dimension array into a list, which is neccessary when exporting the size distribution

    :param array1: the target array
    :type array1: narray
    :return: the transformed list
    
    """

    output=[]
    for i in range(len(array1)):
        output.append(float(array1.item(i)))

    return output

def output_list(string_list, file):
    """
    Export a list of string as one file

    :param string_list: list of string
    :type string_list: list
    :param file: output file name
    :type file: string

    """

    tmp1=open(file,'w',encoding='utf-8')
    tmp1.write('\n'.join(string_list))
    tmp1.close()

def print_list(string_list):
    """
    Display the content of a string list

    :param string_lsit: list of string
    :type string_list: list

    """

    for i in range(len(string_list)):
        print(string_list[i])