# This is the master file

import numpy as np
from scipy.optimize import *
import sys

# Define conversion functions from AES to SI

def convert_inT_ci(inT_ci):
    T_ci=(inT_ci+459.67)*(5/9)
    return T_ci
def convert_inT_hi(inT_hi):
    T_hi=(inT_hi+459.67)*(5/9)
    return T_hi
def convert_T_cout(T_cout):
    T_cout=(T_cout+459.67)*(5/9)
    return T_cout
def convert_T_hout(T_hout):
    T_hout=(T_hout+459.67)*(5/9)
    return T_hout
def convert_inM_c(inM_c):
    M_c=inM_c*.45359
    return M_c
def convert_inM_h(inM_h):
    M_h=inM_h*.45359
    return M_h
def convert_inU(inU):
    U=inU*1055/(.092903*(9/5))
    return U



# Load Thermo data into array

'''
# We can use the following code if we want to load the values into the arrays from the .dat file

Thermo_data = np.loadtxt('Thermophysical_Properties.dat',delimiter=',')

Water_array=np.array(Thermo_data[0,:])
R134a_array=np.array(Thermo_data[1,:])
Ethanol_array=np.array(Thermo_data[2,:])
trimethylpentane_array=np.array(Thermo_data[3,:])
'''

Water_array=np.array([18.01528,273.15,373.15,2.7637E+05,-2.0901E+03,8.1250E+00,-1.4116E-02,9.3701E-06])
R134a_array=np.array([102.03089,172.00,247.08,6.5108E+05,-9.5057E+03,6.2835E+01,-1.8264E-01,2.0031E-04])
Ethanol_array=np.array([46.06844,159.05,351.44,1.0264E+05,-1.3963E+02,-3.0341E-02,2.0386E-03,0.0000E+00])
trimethylpentane_array=np.array([114.22852,165.777,372.388,9.5275E+04,6.9670E+02,-1.3765E+00,2.1734E-03,0.0000E+00])



# Choose fluid type for the cold stream and define all fluid properties

correct=False
while not (correct):
    cfluid_type=input('Enter the fluid type for the cold stream (Water, R134a, Ethanol, 2,2,4-trimethylpentane): ')
    if cfluid_type == "Water":
        cmolecular_weight=Water_array[0] #g/mol
        cmelting_point=Water_array[1] #K
        cboiling_point=Water_array[2] #K
        Ac=Water_array[3]
        Bc=Water_array[4]
        Cc=Water_array[5]
        Dc=Water_array[6]
        Ec=Water_array[7]
        correct=True
    elif cfluid_type == "R134a":
        cmolecular_weight=R134a_array[0]
        cmelting_point=R134a_array[1] #K
        cboiling_point=R134a_array[2] #K
        Ac=R134a_array[3]
        Bc=R134a_array[4]
        Cc=R134a_array[5]
        Dc=R134a_array[6]
        Ec=R134a_array[7]
        correct=True
    elif cfluid_type == "Ethanol":
        cmolecular_weight=Ethanol_array[0]
        cmelting_point=Ethanol_array[1] #K
        cboiling_point=Ethanol_array[2] #K
        Ac=Ethanol_array[3]
        Bc=Ethanol_array[4]
        Cc=Ethanol_array[5]
        Dc=Ethanol_array[6]
        Ec=Ethanol_array[7]
        correct=True
    elif cfluid_type == "2,2,4-trimethylpentane":
        cmolecular_weight=trimethylpentane_array[0]
        cmelting_point=trimethylpentane_array[1] #K
        cboiling_point=trimethylpentane_array[2] #K
        Ac=trimethylpentane_array[3]
        Bc=trimethylpentane_array[4]
        Cc=trimethylpentane_array[5]
        Dc=trimethylpentane_array[6]
        Ec=trimethylpentane_array[7]
        correct=True
    else:
        print("You've entered an invalid fluid type. This program is terminating.")
        sys.exit()



# Choose fluid type for the hot stream and define all fluid properties

correct=False
while not (correct):
    hfluid_type=input('Enter the fluid type for the hot stream (Water, R134a, Ethanol, 2,2,4-trimethylpentane): ')
    if hfluid_type == "Water":
        hmolecular_weight=Water_array[0] #g/mol
        hmelting_point=Water_array[1] #K
        hboiling_point=Water_array[2] #K
        Ah=Water_array[3]
        Bh=Water_array[4]
        Ch=Water_array[5]
        Dh=Water_array[6]
        Eh=Water_array[7]
        correct=True
    elif hfluid_type == "R134a":
        hmolecular_weight=R134a_array[0] #g/mol
        hmelting_point=R134a_array[1] #K
        hboiling_point=R134a_array[2] #K
        Ah=R134a_array[3]
        Bh=R134a_array[4]
        Ch=R134a_array[5]
        Dh=R134a_array[6]
        Eh=R134a_array[7]
        correct=True
    elif hfluid_type == "Ethanol":
        hmolecular_weight=Ethanol_array[0] #g/mol
        hmelting_point=Ethanol_array[1] #K
        hboiling_point=Ethanol_array[2] #K
        Ah=Ethanol_array[3]
        Bh=Ethanol_array[4]
        Ch=Ethanol_array[5]
        Dh=Ethanol_array[6]
        Eh=Ethanol_array[7]
        correct=True
    elif hfluid_type == "2,2,4-trimethylpentane":
        hmolecular_weight=trimethylpentane_array[0] #g/mol
        hmelting_point=trimethylpentane_array[1] #K
        hboiling_point=trimethylpentane_array[2] #K
        Ah=trimethylpentane_array[3]
        Bh=trimethylpentane_array[4]
        Ch=trimethylpentane_array[5]
        Dh=trimethylpentane_array[6]
        Eh=trimethylpentane_array[7]
        correct=True
    else:
        print("You've entered an invalid fluid type. This program is terminating.")
        sys.exit()



# Choose AES or SI units and collect all user input

correct=False
while not (correct):
    units=input("Are you entering units in AES or SI?: ")
    if units == "AES":
        inT_ci=float(input("What is the inlet temperature of the cold stream? (F): "))
        inT_hi=float(input("What is the inlet temperature of the hot stream? (F): "))
        inM_c=float(input("What is the mass flow rate of the cold stream? (lbm/sec): "))
        if inM_c < 0:
            print("You entered a negative mass flow rate, which is not valid. This program is terminating.")
            sys.exit()
        inM_h=float(input("What is the mass flow rate of the hot stream? (lbm/sec): "))
        if inM_h < 0:
            print("You entered a negative mass flow rate, which is not valid. This program is terminating.")
            sys.exit()
        inU=float(input("What is the overall heat transfer coefficient? (BTU/s*ft^2*R): "))
        if inU < 0:
            print("You entered a negative overall heat transfer coefficient, which is not valid. This program is terminating.")
            sys.exit()
        output_flow=input("Are you providing the cold or hot outlet temperature? (Type 'cold' or 'hot'): ")
        if output_flow == "cold":
            T_cout=float(input("What is the cold " + cfluid_type + " outlet temperature? (F)"))
        elif output_flow == "hot":
            T_hout=float(input("What is the hot " + hfluid_type + " outlet temperature? (F)"))
        else:
            print("You didn't enter a valid outlet stream. This program is terminating.")
            sys.exit()
        correct=True
    elif units == "SI":
        inT_ci=float(input("What is the inlet temperature of the cold stream? (K): "))
        inT_hi=float(input("What is the inlet temperature of the hot stream? (K): "))
        inM_c=float(input("What is the mass flow rate of the cold stream? (kg/sec): "))
        if inM_c < 0:
            print("You entered a negative mass flow rate, which is not valid. This program is terminating.")
            sys.exit()
        inM_h=float(input("What is the mass flow rate of the hot stream? (kg/sec): "))
        if inM_h < 0:
            print("You entered a negative mass flow rate, which is not valid. This program is terminating.")
            sys.exit()
        inU=float(input("What is the overall heat transfer coefficient?: (J/s*m^2*K)"))
        if inU < 0:
            print("You entered a negative overall heat transfer coefficient, which is not valid. This program is terminating.")
            sys.exit()
        output_flow=input("Are you providing the cold or hot outlet temperature? (Type 'cold' or 'hot'): ")
        if output_flow == "cold":
            T_cout=float(input("What is the cold " + cfluid_type + " outlet temperature? (K): "))
        elif output_flow == "hot":
            T_hout=float(input("What is the hot " + hfluid_type + " outlet temperature? (K): "))
        else:
            print("You didn't enter a valid outlet stream. This program is terminating.")
            sys.exit()
        correct=True



# Convert from AES to SI if necessary

if units == "AES":
    T_ci=convert_inT_ci(inT_ci)
    print(T_ci)
    T_hi=convert_inT_hi(inT_hi)
    M_c=convert_inM_c(inM_c)
    M_h=convert_inM_h(inM_h)
    U=convert_inU(inU)
    if output_flow == "cold":
        T_cout=convert_T_cout(T_cout)
    elif output_flow == "hot":
        T_hout=convert_T_hout(T_hout)
else:
    T_ci=inT_ci
    T_hi=inT_hi
    M_c=inM_c
    M_h=inM_h
    U=inU



# Confirm that the given temperatures are within the liquid phase range
   
if T_ci < cmelting_point:
    print("Your cold inlet temperature is below the melting point of " + cfluid_type + \
    ". This program is terminating. Please restart and enter a higher cold inlet temperature.")
    sys.exit()
elif T_ci > cboiling_point:
    print("Your cold inlet temperature is above the boiling point of " + cfluid_type + \
    ". This program is terminating. Please restart and enter a lower cold inlet temperature.")
    sys.exit()
if T_hi < hmelting_point:
    print("Your hot inlet temperature is below the melting point of " + hfluid_type + \
    ". This program is terminating. Please restart and enter a higher hot inlet temperature.")
    sys.exit()
elif T_hi > hboiling_point:
    print("Your hot inlet temperature is above the boiling point of " + hfluid_type + \
    ". This program is terminating. Please restart and enter a lower hot inlet temperature.")
    sys.exit()
if output_flow == "cold":
    if T_cout < cmelting_point:
        print("Your cold outlet temperature is below the melting point of " + cfluid_type + \
        ". This program is terminating.")
        sys.exit()
    elif T_cout > cboiling_point:
        print("Your cold outlet temperature is above the boiling point of " + cfluid_type + \
        ". This program is terminating.")
        sys.exit()
if output_flow == "hot":
    if T_hout < hmelting_point:
        print("Your hot outlet temperature is below the melting point of " + hfluid_type + \
        ". This program is terminating.")
        sys.exit()
    elif T_hout > hboiling_point:
        print("Your hot outlet temperature is above the boiling point of " + hfluid_type + \
        ". This program is terminating.")
        sys.exit()



# Calculate the unspecified outlet temperature (T_hout or T_cout)

if output_flow == "cold":
    def Q(T_hout):
        T_h=(T_hi+T_hout)/2
        T_c=(T_ci+T_cout)/2
        residual=(M_h*(Ac+Bc*T_h+Cc*T_h**2+Dc*T_h**3+Ec*T_h**4)*(T_hi-T_hout)) - (M_c*(Ac+Bc*T_c+Cc*T_c**2+Dc*T_c**3+Ec*T_c**4)*(T_cout-T_ci))
        return residual
    z = fsolve(Q,200)
    T_hout=z[0]
    print ("The temperature of the hot outlet is: " + str(T_hout) + "K")
elif output_flow == "hot":
    def Q(T_cout):
        T_h=(T_hi+T_hout)/2
        T_c=(T_ci+T_cout)/2
        residual=(M_h*(Ah+Bh*T_h+Ch*T_h**2+Dh*T_h**3+Eh*T_h**4)*(T_hi-T_hout)) - (M_c*(Ah+Bh*T_c+Ch*T_c**2+Dh*T_c**3+Eh*T_c**4)*(T_cout-T_ci))
        return residual
    z = fsolve(Q,200)
    T_cout = z[0]
    print ("The temperature of the cold outlet is: " + str(T_cout) + "K")



# Use the provided equations to solve for the other unknowns
# (I've temporarily added print statements to make sure each value is correct)

delta_T2 = T_hout-T_ci
print("delta_T2 is: " + str(delta_T2))
delta_T1 = T_hi-T_cout
print("delta_T1 is: " + str(delta_T1))

delta_Tlm = (delta_T2-delta_T1)/(np.log(delta_T2/delta_T1))
print("delta_Tlm is: " + str(delta_Tlm))

R = (T_hi-T_hout)/(T_cout-T_ci)
print("R is: " + str(R))
P = (T_cout-T_ci)/(T_hi-T_ci)
print("P is: " + str(P))

F = (np.sqrt(R**2+1)/(R-1)) * ((np.log((1-P)/(1-P*R))) / (np.log((2-P*(R+1-np.sqrt(R**2+1)))/(2-P*(R+1+np.sqrt(R**2+1))))))
print("F is: " + str(F))

T_h=(T_hi+T_hout)/2
print("T_h is: " + str(T_h))
T_c=(T_ci+T_cout)/2
print("T_c is: " + str(T_c))


# Divide the heat capacities by 1000 to convert from kmol to mol

C_ph = (Ah+Bh*T_h+Ch*T_h**2+Dh*T_h**3+Eh*T_h**4)/1000
print("C_ph is: " + str(C_ph))
C_pc = (Ac+Bc*T_c+Cc*T_c**2+Dc*T_c**3+Ec*T_c**4)/1000
print("C_pc is: " + str(C_pc))


# Divide the molecular weight by 1000 to convert from g/mol to kg/mol

qh = M_h*(C_ph/(hmolecular_weight/1000))*(T_hi-T_hout)
print("qh is: {:.3e}".format(qh) + " W")
qc = M_c*(C_pc/(cmolecular_weight/1000))*(T_cout-T_ci)
print("qc is: {:.3e}".format(qc) + " W")

A = qh/(F*U*delta_Tlm)
print("The ideal surface area is: {:.2f}".format(A) + " m^2")
A = qc/(F*U*delta_Tlm)
print("The ideal surface area is: {:.2f}".format(A) + " m^2")

Cost = A*1000
print("The cost of the custom heat exchanger is: $" + str(Cost))

qt = F*U*A*delta_Tlm
print("qt is: {:.3e}".format(qt) + " W")
