import numpy as np
import math

def func(x1,x2):            
    f1 = (1+(x1+x2-2*(x1**2)-(x2**2)+x1*x2))    # Fitness Function
    func_max = 1/f1
    return func_max

def inv_func(f):
    return 1-(1/f)      # Inverse of the Fitness Function
    
x1_max = x2_max = 0.5    # Maximum Value of the Varialbles x1, x2
x1_min = x2_min = 0.0    # Minimum Value of the Varialbles x1 ,x2

Pc = 1        # Crossover Probability    
Gen=150       # Number Of Generations
n=20            # String Size of each variable
N=6             # Population Size of Solution Matrix

avg_fit_eval = np.zeros([Gen,2])
max_fit_val = np.zeros([Gen,3])
x_opt = np.zeros([Gen,3])
inv_fit_val = np.zeros([Gen,2])

S = np.zeros([N,2*n])         # Solution Matrix
dec = np.zeros([n,1])          # Decoding Vector
print(dec)

x = np.zeros([N,2])
x_norm = np.zeros([N,2])    # Normalized x
fit_eval = np.zeros([N,1])          # Fitness Evaluation

cso=0                       # Crossover Count
for i in range (0,n-1):
    dec[i]=pow(2,n-i-1)      # Binary String Decoding Vector
print(dec)

# Initial Fitness Values
for i in range (0,N):
    x1=np.random.randint(2, size=n)
    x2=np.random.randint(2, size=n)
    s=np.concatenate((x1,x2))
    S[i,:]=s
    for j in range (0,n-1):
        x[i,0]=x[i,0]+dec[j]*x1[j]
        x[i,1]=x[i,1]+dec[j]*x2[j]
    x_norm[i,0]=x1_min+((x1_max-x1_min)/(pow(2,n)-1))*x[i,0]
    x_norm[i,1]=x2_min+((x2_max-x2_min)/(pow(2,n)-1))*x[i,1]
    fit_eval [i,0]=func(x_norm[i,0],x_norm[i,1])       # Initial Fitness Evaluation
# print(fit_eval )

print("Gen \t optimam_x1 \t\t optimam_x2 \t\t Avg_Fitness \t\t Inv_Fitness")

for g in range (0,Gen):
    sum1=0
    rws_p = np.zeros([N,1])       # Roulette Wheel Selection Probability
    cuv_fit_p = np.zeros([N,1])      # Cumilative Fitness Probability

    for i in range (0,N): 
        sum1=sum1+fit_eval[i,0]
    for i in range (0,N):
        rws_p[i,0]=(fit_eval[i,0]/sum1)
    for i in range (0,N):
        for j in range (0,i+1):
            cuv_fit_p[i,0]=cuv_fit_p[i,0]+rws_p[j,0]
            
    mat_pool = np.zeros([N,2*n])      # Mating Pool

    for i in range (0,N):
        rp = np.random.rand()
        for j in range(0,N):
            if(cuv_fit_p[j,0]>rp):
                mat_pool[i,:]=S[j,:]
                break

    # Using Two Point CrossOver

    ch1 = np.zeros([1,2*n])     # Child Solution 1
    ch2 = np.zeros([1,2*n])     # Child Solution 2

    for i in range (0,int(N/2)):
        Rc = np.random.rand(1)
        if (Rc<Pc):
            rc1=np.random.randint(1,2*n-1)
            rc2=np.random.randint(1,2*n-1)
            rs1=np.random.randint(0,N-1)
            rs2=np.random.randint(0,N-1)

            ch1[0,:]=mat_pool[rs1,:]
            ch2[0,:]=mat_pool[rs2,:]
            if (rc1>rc2):
                ch1[0,rc2:rc1]=mat_pool[rs2,rc2:rc1]
                ch2[0,rc2:rc1]=mat_pool[rs1,rc2:rc1]
            if (rc2>rc1):
                ch1[0,rc1:rc2]=mat_pool[rs2,rc1:rc2]
                ch2[0,rc1:rc2]=mat_pool[rs1,rc1:rc2]
            mat_pool[rs1,:]=ch1[0,:]
            mat_pool[rs2,:]=ch2[0,:]
            cso=cso+1

    # Decoding the fitness Values 
    avg_fit_eval=0
    opt_x1=0
    opt_x2=0
    avg_inv = 0
    fit_eval_new=np.zeros([N,1])       # Updated Fitness Evaluation
    x_norm_new = np.zeros([N,2]) # Updated Normalized x
    x_n = np.zeros([N,2])          # Norm of x      

    for i in range (0,N):
        for j in range (0,2*n):
            if(j<n):
               x_n[i,0]=x_n[i,0]+mat_pool[i,j]*dec[j]
            else:
                x_n[i,1]=x_n[i,1]+mat_pool[i,j]*dec[j-n]
        x_norm_new[i,0]=x1_min+((x1_max-x1_min)/(pow(2,n)-1))*x_n[i,0]
        x_norm_new[i,1]=x2_min+((x2_max-x2_min)/(pow(2,n)-1))*x_n[i,1]
        fit_eval_new[i,0]=func(x_norm_new[i,0],x_norm_new[i,1])
        avg_fit_eval=avg_fit_eval+(1/N)*(fit_eval_new[i,0])
        opt_x1=opt_x1+(1/N)*(x_norm_new[i,0])
        opt_x2=opt_x2+(1/N)*(x_norm_new[i,1])
        avg_inv = avg_inv + (1/N)*(inv_func(fit_eval_new[i,0]))
    max_fit_val[g,0]=g+1
    max_fit_val[g,1]=max(fit_eval_new)
    max_fit_val[g,2]=min(fit_eval_new)
    x_opt[g,0]=g+1
    x_opt[g,1]=opt_x1
    x_opt[g,2]=opt_x2
    inv_fit_val[g,0]=g+1
    inv_fit_val[g,1]=math.fabs(avg_inv)
    fit_eval =fit_eval_new
    S=mat_pool
    print(g+1,"\t",opt_x1,"\t",opt_x2,"\t",avg_fit_eval,"\t",math.fabs(avg_inv))


print("Average Fitness Value : ",avg_fit_eval)
print("\nNumber of Cross Over Operations : ",cso)
print("\n\nNo. of Gen\tMax Fitness\tMin Fitness Values \n ",max_fit_val)



















