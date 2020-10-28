"""
Created on Oct 2020

@author: murina
"""

import numpy as np



def axis_for_graph(_l0, l0_ , _l1, l1_, _l2, l2_, step, seats, persons, l, m):
    
    x = []

    y = []

    z = []

    ex = []

    ey = []

    ez = []

    e_total = []
    
    bit_sum = []

    for lambda_0 in np.arange( _l0, l0_, step):
        for lambda_1 in np.arange( _l1, l1_, step):
            for lambda_2 in np.arange( _l2, l2_, step):

                a, b, c, d, e = tuples_for_graph(seats, persons, lambda_0, lambda_1, lambda_2, l, m)
                x.append(lambda_0)
                y.append(lambda_1)
                z.append(lambda_2)
                e_total.append(b)
                ex.append(c)
                ey.append(d)
                ez.append(e)
                bit_sum.append(a)

    return bit_sum,x,y,z,ex,ey,ez, e_total



def tuples_for_graph(seats, persons, lambda_0, lambda_1, lambda_2, l, m):
    
    
    import dimod   
    
    from collections import defaultdict 

    h = defaultdict(float)

    J = defaultdict(float)
    
    e_offset = lambda_0 * persons**2
       
    for _ in range(seats):
        
        h[ _ ] = lambda_0 * ( 1 - 2*persons )   
    
    for _ in range(seats - 1):
        for __ in range( _ + 1, seats):
            J[ _ , __ ] =  lambda_0*2    
        
    
     
        
    for _ in range(seats-1):
        J[ _ , _ + 1] += lambda_1 * m
        
    J[0, seats - 1] += lambda_1 * m
    
    
    
    
    quarter_period = int(seats/4)
    
    for _ in range(quarter_period):
        J[ _ , 3*quarter_period -1 - _ ] = lambda_2 * l
        
    for _ in range(quarter_period, 2*quarter_period):
        J[ _ , 5*quarter_period-1 - _ ] = lambda_2 * l
    
       
    
    bqm = dimod.BinaryQuadraticModel(h,J,e_offset,'BINARY')   #dimod.BINARY)
    
    results = dimod.ExactSolver().sample(bqm)
    
    
    
    smpl = results.first.sample
    

  
    
    e_lambda0 = 0.
    
    for _ in range(len(smpl)):
        e_lambda0 += smpl[_]
    
    e_lambda0 += -persons 
    
    energy_e_lambda0 = lambda_0 * e_lambda0**2 #lambda_0 * ( e_lambda0**2 - persons**2 )
    
    
    
    
    e_lambda1 = smpl[0]*smpl[len(smpl)-1]
    
    
    for _ in range(len(smpl)-1):
        e_lambda1 += smpl[_]*smpl[_+1]
        
    
    energy_e_lambda1 = lambda_1*m*e_lambda1 
    
   
    
    e_lambda2 = 0
    
    for _ in range(quarter_period):
        e_lambda2 += smpl[ _ ] * smpl[ 3*quarter_period -1 - _ ]  
        
    for _ in range(quarter_period, 2*quarter_period):
        e_lambda2 += smpl[ _ ] * smpl[ 5*quarter_period-1 - _ ]
        
    energy_e_lambda2 = lambda_2*l*e_lambda2 


    total = energy_e_lambda0 + energy_e_lambda1 + energy_e_lambda2 

    bit_sum = 0 

        
    for key, values in smpl.items():
        bit_sum += values
 
        
    
    
    return  bit_sum, total, energy_e_lambda0, energy_e_lambda1, energy_e_lambda2


def ploting( x , y, z, ex, ey, ez, e_total):

    import matplotlib.pyplot as plt
    
    import numpy as np
    
    
    axis = np.arange(len(x))
    
    fig, ax = plt.subplots(3,2)
    
    fig.set_figheight(15)
    fig.set_figwidth(15)
    
    ax[0,0].set(ylabel='lambda 0')
    ax[1,0].set(ylabel='lambda 1')
    ax[2,0].set(ylabel='lambda 2')
    ax[0,1].set(ylabel='energy 0')
    ax[1,1].set(ylabel='energy 1')
    ax[2,1].set(ylabel='energy 2')
    
    ax[0,0].plot(axis, x)
    ax[1,0].plot(axis, y)
    ax[2,0].plot(axis, z)
    ax[0,1].plot(axis, ex)
    ax[1,1].plot(axis, ey)
    ax[2,1].plot(axis, ez)
    
    fig2, ax2 = plt.subplots(1,1)
    
    fig2.set_figwidth(20)
    
    ax2.set(title='total energy = energy 0 + energy 1 + energy 2')
    ax2.scatter(axis, e_total)



