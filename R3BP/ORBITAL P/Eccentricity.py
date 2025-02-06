import matplotlib.pyplot as plt
import numpy as np

#ORBPAR= eccentricity,semi-majoraxis,long.ascendingnode,inclination,argumentofperiapsis 
#5graphstotal

num_observationsH=[3,6,9,12,15,18,21,24,27,30]
helio_uncertainty = [0.001570363308,0.001571844221,0.001571844221,0.001574826287,0.00157294234,0.00157813036,0.00157813036,0.001585783949,0.001589206577,0.00159247595]  

num_observationsK= [3,6,9,12,15,18,21,24,27,30]
kallisto_uncertainty = [0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01]  

num_observationsA=[3,6,9,12,15,18,21,24,27,30]
aneas_uncertainty = [0.02081167404,0.00938346858,0.00938346858,0.01353782977,0.02338115986,0.02110302894,0.02110302894,4.27E-05,3.01E-05,2.95E-05]  

num_observationsM=[]
marion_uncertainty=[]

num_observationsD=[]
deborah_uncertainty=[]

num_observationsT=[]
themis_uncertainty=[]

num_observationsE=[]
eurykleia_uncertainty=[]



plt.figure(figsize=(8, 6))
plt.title('Fractional Uncertainty vs. Number of Observations')
plt.plot(num_observationsH, helio_uncertainty, label='Asteroid 1', marker='o')
plt.plot(num_observationsK, kallisto_uncertainty, label='Asteroid 2', marker='s')
plt.plot(num_observationsA, aneas_uncertainty, label='Asteroid 3', marker='^')
plt.plot(num_observationsM, marion_uncertainty, label='Asteroid 4', marker='o')
plt.plot(num_observationsT, themis_uncertainty, label='Asteroid 5', marker='s')
plt.plot(num_observationsD, deborah_uncertainty, label='Asteroid 6', marker='^')
plt.plot(num_observationsE, eurykleia_uncertainty, label='Asteroid 7', marker='o')
plt.xlabel('Number of Observations')
plt.ylabel('Fractional Uncertainty')
plt.grid(True)
plt.legend()
plt.show()






