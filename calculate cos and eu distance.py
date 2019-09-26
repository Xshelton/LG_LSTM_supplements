import numpy as np
cos= np.load('cos_list.npy') 
eu=np.load('eu_list.npy')
print(cos)

print(len(cos))
print(len(eu))

print('cos average',np.mean(cos))
print('eu distance average',np.mean(eu))

