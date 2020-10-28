"""
Created on Oct 2020

@author: murina
"""

def run_jobs(seats, persons, lambda_0, lambda_1, lambda_2, l, m, solver):
    
    
    import dimod   
    
    from dwave.system import LeapHybridSampler
    
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
    
    
    if solver == "EXACT":
    
        results = dimod.ExactSolver().sample(bqm)
    
    if solver == "HYBRID":
    
    
    
        results = LeapHybridSampler().sample(bqm)
    
    smpl = results.first.sample
    
    quarter_period = int(seats/4)


    e_lambda0 = 0.
    
    for _ in range(len(smpl)):
        e_lambda0 += smpl[_]
    
    
    
    e_lambda0 += -persons 
    
    
    energy_e_lambda0 = lambda_0 * e_lambda0**2 
    
    
   
    
    
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

    

    ham = "Output configuration:      >"
    
    oruxo = 0
    
    for _ , __ in smpl.items():
        ham += str(__) + "  "
        oruxo += __
    
        
    print("Problem dimension (# of variables -seats-) = ", seats )
    print("Number of knigths = ", persons )

    print(ham)    
    print("")
    print("")
    
    spam = "            "
    
    egg = "          "
    
    for _ in reversed(range(2*quarter_period, 3*quarter_period)):
        spam += str(smpl[_]) + "     "
        egg += "* * * "
    print(spam) 
    print(egg)
    
    print("       *", quarter_period*"      ", "*  ")
    
    egg = "          "
    
    for _ in range(3*quarter_period, 4*quarter_period):
        print(smpl[_],"     *", quarter_period*"      ", "*  ",smpl[5*quarter_period - _ -1])   
        print("       *", quarter_period*"      ", "*  ")  
        egg += "* * * "
    
    print("       *", quarter_period*"      ", "*  ")
    print(egg)
    
    spam = "           >"
    
    for _ in range(quarter_period):
        spam += str(smpl[_]) + "     "
    print(spam) 
    print("")
    print("Number of sitting knigths = ", oruxo )
    
    
    
    
    print("")
    print("")
    print("Energy for the output configuration")
    print("_______________________________________________________________________________")
    print("Constraint I (all the knights sitting) --------------------------->  =", "%.2f" % energy_e_lambda0)
    print("Constraint II (no two knigths sitting next to each other) -------->  =", "%.2f" % energy_e_lambda1)
    print("Constraint III (no two knigths sitting one in front the other) --->  =", "%.2f" % energy_e_lambda2)
    print("                                                                    -----------")
    print("                                                               Total =", "%.2f" % total)
    print("")    
    print("lambda_0  |  lambda_1  |  lambda_2   ")
    print("  %.2f   " % lambda_0,"|", "  %.2f    " % lambda_1,"|" ,"   %.2f" % lambda_2 )
    
    return  bit_sum