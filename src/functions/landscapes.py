import numpy as np
from matplotlib import pyplot as plt
import math 
from scipy.interpolate import interp1d
from scipy import integrate

# Sorts the list to be increasing in first variable, then decreasing in the second variable
def pl_sort(A):
    A = [(x[0],x[1]) for x in A]
    A.sort()
    a=np.array(A)
    sep = np.split(a, np.where(np.diff(a[:,0]))[0]+1)
    flip = [np.flip(x, axis=0) for x in sep]
    S = np.concatenate(flip)
    lst = [(x[0],x[1]) for x in S]
    return lst

# Formerly pl
def persistence_landscape(barcode): #input: a list of 2-tuples. be careful if np array!
    a = pl_sort(barcode)
    PL=[]
    while len(a)>0:
        p=0
        (b,d) = a.pop(0)
        L = []
        L.extend([(-np.inf,0),(b,0),((b+d)/2,(d-b)/2)])
        while L[-1] != (np.inf,0):
            ds = [x[1] for x in a[p:]]
            if len(ds)==0:
                L.extend([(d,0),(np.inf,0)])
            elif d >= max(ds):
                L.extend([(d,0),(np.inf,0)])
            else:
                for j in range(p,len(a)):
                    if a[j][1]>d:
                        (B,D) = a[j]
                        break
                        break
                    else:
                        continue 
                J = a.index((B,D)) 
                a.pop(J)
                p = J  
                if B>d:
                    L.append((d,0))
                if B>=d:
                    L.append((B,0))
                else:
                    L.append(((B+d)/2,(d-B)/2))
                    a.append((B,d))
                    a = pl_sort(a)
                    p = a.index((B,d)) +1 #unsure whether this index is right 
                L.append(((B+D)/2,(D-B)/2))
                (b,d) = (B,D)
        PL.append(L)
    return PL

# formerly pl_plot
def plot_pl(pl, xl = 2, xu = 12, yl = 0, yu = 4, title = ""):
    for item in pl:
        x = [i[0] for i in item]
        y = [i[1] for i in item]
        plt.scatter(x,y,s=1)
        plt.plot(x,y)
        plt.xlim(xl,xu)
        plt.ylim(yl,yu)
        plt.title(title)
        plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

# make list of all x coords in each PL degree for two lists of PL diagrams, without repetition 
# formerly make_X
def common_x(a,b): 
    x_a = [[i[0] for i in item[1:-1]] for item in a] #remove inf for sorting but add back later
    x_b = [[i[0] for i in item[1:-1]] for item in b]
    n = min(len(a),len(b))
    X = [x+y for x,y in zip(x_a[:n],x_b[:n])]
    if len(a)>len(b):
        X.extend(x_a[n:])
    else:
        X.extend(x_b[n:])
    X = [sorted(list(set(x))) for x in X]
    for i in range(len(X)):
        X[i].insert(0,-math.inf) #add inf and -inf back 
        X[i].append(math.inf)
    return X

# turn the persistence landscape a into an interpolated function
# formerly pl2func
def pl_to_function(pl):
    x = [[k[0] for k in item] for item in pl] 
    y = [[k[1] for k in item] for item in pl]
    L = []
    for i in range(len(pl)):
        l = interp1d(x[i],y[i],kind = 'linear')
        L.append(l)
    return L


# formerly avg_2pl
def average_pair(pl1,pl2): #take as input 2 persistence landscapes
    a=pl1
    b=pl2
    X = common_x(a,b)
    n = min(len(a),len(b))
    N = max(len(a),len(b))
    
    L_a = pl_to_function(a)
    L_b = pl_to_function(b)
  
    Y_a = [] #values a takes on all of X- needed the interpolated functions to find these
    for i in range(len(a)):
        y = [L_a[i](x) for x in X[i]]
        Y_a.append(y)
    if len(a)<len(b):
        for i in range(n,N):
            Y_a.append([0 for x in X[i]])

    Y_b = [] #values b takes on all of X
    for i in range(len(b)):
        y = [L_b[i](x) for x in X[i]]
        Y_b.append(y)
    if len(b)<len(a):
        for i in range(n,N):
            Y_b.append([0 for x in X[i]])
        
    X = [np.array(x) for x in X] #average them - note we only go up to the highest k that they both have pls for
    Y_a = [np.reshape(item,len(item)) for item in Y_a] #reshape into numpy arrays 
    Y_b = [np.reshape(item,len(item)) for item in Y_b]
    Y_avg = [np.mean([Y_a[i],Y_b[i]], axis = 0) for i in range(N)] #these are just the Y coords, need to get it back into a pl 
    
    avg_pl = []
    for i in range(N):
        pair = list(zip(X[i],Y_avg[i]))
        avg_pl.append(pair)
    
    return avg_pl   

# formerly avg_pl
def average_list(lst_pl): #input a list of persistence landscapes
    old_l = lst_pl
    while len(old_l)>1:
        n = len(old_l)
        new_l = []
        if n % 2 == 0:
            for i in range(int(n/2)):
                new_l.append(average_pair(old_l[2*i],old_l[2*i+1]))
        else:
            for i in range(int((n-1)/2)):
                new_l.append(average_pair(old_l[2*i],old_l[2*i+1]))
            new_l.insert(int((n-1)/2),old_l[-1])
        old_l = new_l
    return old_l[0]    

def classic_integration(pl,k): #integrates kth func in the persistence landscape
    x_vals = [[i[0] for i in item] for item in pl]
    L_a = pl_to_function(pl)
    result = 0
    for i in range(len(x_vals[k])-1):
        (x,y) = (x_vals[k][i],x_vals[k][i+1])
        integr = integrate.quad(lambda x: L_a[k](x),x,y)[0]
        result += integr
    return result

# trapezium integration
# formerly tp_int_pl
def trapezium_integrate(pl,k): #integrates kth func in the persistence landscape
    x_vals = [[i[0] for i in item][1:-1] for item in pl] #get rid of infinities
    L_a = pl_to_function(pl)
    result = 0
    for i in range(len(x_vals[k])-1): #need to get rid of the infinities
        (x,y) = (x_vals[k][i],x_vals[k][i+1])
        trap = (y-x)*(L_a[k]((1/2)*(y+x)))
        result += trap
    return result

def add_pl(pl1,pl2): #only adds in the dimensions that both have persistence landscapes
        a=pl1
        b=pl2
        X = common_x(a,b)
        L_a = pl_to_function(a)
        L_b = pl_to_function(b)

        Y_a = [] #values a takes on all of X- needed the interpolated functions to find these
        for i in range(len(a)):
            y = [L_a[i](x) for x in X[i]]
            Y_a.append(y)

        Y_b = [] #values b takes on all of X
        for i in range(len(b)):
            y = [L_b[i](x) for x in X[i]]
            Y_b.append(y)

        X = [np.array(x) for x in X] #average them - note we only go up to the highest k that they both have pls for
        Y_a = [np.reshape(item,len(item)) for item in Y_a] #reshape into numpy arrays 
        Y_b = [np.reshape(item,len(item)) for item in Y_b]
        n = min(len(a), len(b)) #only average over degrees they have in common 
        Y_ab = [np.add(Y_a[i],Y_b[i]) for i in range(n)] #these are just the Y coords, need to get it back into a pl 

        sm = []
        for i in range(n):
            pair = list(zip(X[i],Y_ab[i]))
            sm.append(pair)
        
        return sm   


# formerly n_avg
def average_n(lst_pl):
    sm = lst_pl[0]
    for i in range(len(lst_pl)-1):
        sm = add_pl(sm,lst_pl[i+1])
    n = len(lst_pl)
    avg = [[(i[0],i[1]/n) for i in item] for item in sm]
    return avg

def L_inf(p1,p2): #given two persistence landscapes returns a distance
    a = p1
    b = p2
    n = min(len(a),len(b))

    X = common_x(a,b) #make list of points X union of x coords of critical pts in both 
    
    L_a = pl_to_function(a)
    L_b = pl_to_function(b)
        
    Y_a = [] #Y values a takes on all of X, but only up to min length 
    for i in range(n):
        y = [L_a[i](x) for x in X[i]]
        Y_a.append(y)

    Y_b = [] #values b takes on all of X, but only up to min length 
    for i in range(n):
        y = [L_b[i](x) for x in X[i]]
        Y_b.append(y)
        
    X = [np.array(x) for x in X]
    Y_a = [np.reshape(item,len(item)) for item in Y_a]
    Y_b = [np.reshape(item,len(item)) for item in Y_b]
    
    Y_diff = [item1-item2 for item1,item2 in zip(Y_a,Y_b)] #this is a list of lists of scalars, their difference function NOT NECC A PL 
    norms = [[np.abs(i) for i in item] for item in Y_diff] 
    distance = max([max(i) for i in norms])
    return distance

def Lp(pl1,pl2,p):
    a=pl1
    b=pl2
    n = min(len(a),len(b))
    X = common_x(a,b)
    L_a = pl_to_function(a)
    L_b = pl_to_function(b)
    
    Y_a = [] #Y values a takes on all of X
    for i in range(len(a)):
        y = [L_a[i](x) for x in X[i]]
        Y_a.append(y)

    Y_b = [] #values b takes on all of X
    for i in range(len(b)):
        y = [L_b[i](x) for x in X[i]]
        Y_b.append(y)
        
    X = [np.array(x) for x in X]
    Y_a = [np.reshape(item,len(item)) for item in Y_a]
    Y_b = [np.reshape(item,len(item)) for item in Y_b]
    Y_diff_abs = [np.absolute(item1-item2)**p for item1,item2 in zip(Y_a,Y_b)]
    diff_pl = [] #this is our new pl
    for i in range(n):
        pair = list(zip(X[i],Y_diff_abs[i]))
        diff_pl.append(pair)
    diff_func = pl_to_function(diff_pl) #turn it into a function (actually it is a list of functions) dont need this yet
    
    sum_int = 0
    for i in range(len(diff_pl)):
        j = trapezium_integrate(diff_pl,i) #use my piecewise trapezium integration function
        sum_int += j
        
    final = (sum_int)**(1/p)
    return final   


# formerly dm_lp
def lp_distance_matrix(pls,p):  
    all_distances = np.zeros((len(pls), len(pls)))
    for i in range(0, len(pls)):
        for j in range(i, len(pls)):
            all_distances[i,j] = Lp(pls[i], pls[j],p) 
    for i in range(0, len(pls)):
        for j in range(0, i):
            all_distances[i,j] = all_distances[j,i]
    return(all_distances)


# formerly dm_linf
def linf_diatance_matrix(pls):    
    all_distances = np.zeros((len(pls), len(pls)))
    for i in range(0, len(pls)):
        for j in range(i, len(pls)):
            all_distances[i,j] = L_inf(pls[i], pls[j]) 
    for i in range(0, len(pls)):
        for j in range(0, i):
            all_distances[i,j] = all_distances[j,i]
    return(all_distances)