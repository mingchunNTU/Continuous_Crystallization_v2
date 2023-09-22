from data import *
from figure import *
from utility import *
import numpy as np
from scipy.optimize import fsolve

class CSD:
    """
    Module used to represent crystal size distribution
    """

    def __init__(self,size,fraction):
        self.size=size
        self.fraction=fraction
        
class PDF:

    """
    Module used to represent population density function
    """
    def __init__(self,size,density_function):
        self.size=size
        self.density_function=density_function

class crystal:

    def read_setting(self,file):
        tmp1=parameter_read(file)
        self.k=get_variable(tmp1,"shape factor").value[0]
        self.rho=get_variable(tmp1,"crystal density").value[0]
        self.Kr=get_variable(tmp1,"Kr").value[0]
        self.j=get_variable(tmp1,"j").value[0]
        self.i=get_variable(tmp1,"i").value[0]
        self.tau=get_variable(tmp1,"residence time").value[0]
        self.MT=get_variable(tmp1,"slurry concentration").value[0]

    def nucleation(self,G):
        G_used=G*1e-6/60 # um/min to m/s
        B=self.Kr*(self.MT**self.j)*(G_used**self.i)*1000*1e-18*60 # 1/L*s to 1/um^3*min
        return B

    def material_balance(self,G):
        n=self.nucleation(G)/G
        MT=6*self.k*self.rho*n*(G*self.tau)**4
        output=MT-self.MT
        return output
    
    def solve(self):
        G=fsolve(self.material_balance,[5])
        
        self.G=G
        self.B=self.nucleation(G)
        self.n=self.B/G
        self.MT2=6*self.k*self.rho*self.n*(G*self.tau)**4

    def output(self,setting_dir,CSD_mesh_size):
        """
        Export the simulation result (PDF and CSD)
        """

        number=int(15*self.G*self.tau/CSD_mesh_size)+1
        size1=variable("size","$\mu m$",array_to_list(np.linspace(0,15*self.G*self.tau,number)))
        density_function=variable("size","$1/\mu m^4$",[])
        for i in range(len(size1.value)):
            n=float(self.n*np.exp(-size1.value[i]/self.G/self.tau))
            density_function.value.append(n)

        self.product_PDF=PDF(size1,density_function)

        # export the PDF
        PDF_output=[self.product_PDF.size,self.product_PDF.density_function]
        variable_output(PDF_output,setting_dir+"Result/PDF.csv")

        # export the CSD
        mesh_size=CSD_mesh_size
        size,number_fraction=DF_to_NF(self.product_PDF.size,self.product_PDF.density_function,mesh_size)
        CSD_output=[size,number_fraction]
        variable_output(CSD_output,setting_dir+"Result/CSD(number).csv")

        size,volume_fraction=NF_to_VF(size,number_fraction)
        CSD_output=[size,volume_fraction]
        variable_output(CSD_output,setting_dir+"Result/CSD(volume).csv")


        # print the simulation result
        tmp1=[]
        tmp1.append("Growth Rate: "+str(self.G)+" um/min")
        tmp1.append("Nucleation Rate: "+str(self.B)+" 1/(um^3*min)")
        tmp1.append("Residence Time: "+str(self.tau)+" min")
        tmp1.append("G*tau: "+str(self.G*self.tau)+" um")
        tmp1.append("MT2: "+str(self.MT2)+" kg/m3")
        print_list(tmp1)
        output_list(tmp1,setting_dir+"Result/report.txt")
