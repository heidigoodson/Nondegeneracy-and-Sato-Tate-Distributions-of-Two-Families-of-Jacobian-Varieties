"""
This file stores the functions needed to perform computations in the paper [1].
The algorithms are based on [1] Emory and Goodson "Nondegeneracy and Sato-Tate Distributions of Two Families of Jacobian Varieties" arXiv 2024

When using this code in your research, please cite as:
"the software of [1], where [1] is as above."


"""


##########################################################################

#Necessary building blocks
I = identity_matrix(QQ,2)
O = zero_matrix(QQ,2)
J = matrix(QQ,[[0,1],[-1,0]])


##########################################################################

#Functions for the curves y^2=x^n-x where n=2^d+1

#These functions build the gamma matrix
def gammaijx(i,j,d):
    a=5 #generator
    modulus=2^(d+1)
    g = Integer(2^(d-1)) #genus of the curve
    if 2*j-1==Integer(mod(a*(2*i-1),modulus)):
        return(I)
    elif j<=(g/2) and 2*j-1==modulus-Integer(mod(a*(2*i-1),modulus)):
        return (-J)
    elif j>(g/2) and 2*j-1==modulus-Integer(mod(a*(2*i-1),modulus)):
        return (J)
    else:
        return (O)

def gammarowsx(d): #create a list of rows of the block matrix
    g = Integer(2^(d-1)) #genus of the curve
    rows=[]
    for i in [1..g]:
        rowi = []
        for j in [1..g]:
            rowi.append(gammaijx(i,j,d))
        rows.append(rowi)
    return(rows)

def gammamatrix_x(d): #working from list of rows, build the gamma matrix as a block matrix with g block rows and g block columns
    g = Integer(2^(d-1)) #genus of the curve
    if d<2:
        return 'Choose d>1'
    return(block_matrix(g,g,gammarowsx(d)))

def matrixoutputx(d):
    print("This is the gamma matrix for the curve y^2=x^",2^d+1,"-x")
    print("")
    print(gammamatrix_x(d))

#Checking the order of gamma in the component group
def GammaOrderx(d):
    if d<2:
        return 'Choose d>1'
    gamma=gammamatrix_x(d)
    maybeorder=Integer(euler_phi(2^(d+1))/2) #this is what I think is the order of gamma
    g = Integer(2^(d-1)) #genus of the curve
    for n in [1..maybeorder]:
        if gamma^n==identity_matrix(2*g) or gamma^n==-identity_matrix(2*g):
            return n
    return 0

#This function builds the other generator of the component group
#order=2
def Jmatrix_x(d):
    g = Integer(2^(d-1)) #genus of the curve
    if d<2:
        return 'Choose d>1'
    blocks=[]
    for i in [1..g]:
        blocks.append(J)
    return block_diagonal_matrix(blocks)

#This function builds the generator of the component group for generic c

def gammacx(d):
    g = Integer(2^(d-1)) #genus of the curve
    F.<z>=CyclotomicField(4*g)
    Zd=diagonal_matrix([z,z])
    Zblocks=[]
    for i in [0..g-1]:
        Zi=Zd^(2*i-g+1)
        Zblocks.append(Zi)
    return block_diagonal_matrix(Zblocks)


#identity component for the curve y^2=x^n-x, n=2^d+1
def IDx(d):
    deg=2^(d+1)
    if d<2:
        return 'Choose d>1'
    F.<z>=CyclotomicField(deg)
    g=Integer(2^(d-1))
    P=PolynomialRing(F,g/2,'u')
    blocks=[]
    for i in [0..(g/2)-1]:
        blocks.append(diagonal_matrix([P.gens()[i],(P.gens()[i])^(-1)]))
    return block_diagonal_matrix(blocks+blocks)


#This prints the characteristic polynomials for a given d
def charpolyx(d, coeff=1):
    if d<2:
        return 'Choose d>1'
    print("Here are the characteristic polynomials for y^2=x^",2^d+1,"-x")
    print("To save space, it only prints the a_n coefficient, with default a_1")
    print("")
    g = Integer(2^(d-1)) #genus of the curve
    G0=IDx(d)
    gamma=gammamatrix_x(d)
    gammaorder=GammaOrderx(d) #this is the order of gamma!
    Jd=Jmatrix_x(d) #the order of this in ST/ST^0 is always 2
    for k in [0..gammaorder-1]:
        for j in [0..1]:
            f=(G0*gamma^k*Jd^j).characteristic_polynomial()
            nice_f=(f(SR('x'))).expand().collect(x) #this makes the coefficients look a little nicer
            print("k=",k,"j=",j, f[coeff])
            #print("k=",k,"j=",j,nice_f[coeff])
            print("")


#generic c

def charpolycx(d,coeff=1):
    if d<2:
        return 'Choose d>1'
    print("Here are the characteristic polynomials for y^2=x^",2^d+1,"-cx")
    print("To save space, it only prints the a_n coefficient, with default a_1")
    print("")
    g = Integer(2^(d-1)) #genus of the curve
    G0=IDx(d)
    gamma=gammamatrix_x(d)
    gammaorder=GammaOrderx(d) #this is the order of gamma!
    Jd=Jmatrix_x(d) #the order of this in ST/ST^0 is always 2
    gammac=gammacx(d) #the order is 2g
    for i in [0..2*g-1]:
        f=(G0*gamma^0*Jd^0*gammac^i).characteristic_polynomial()
        nice_f=(f(SR('x'))).expand().collect(x) #this makes the coefficients look a little nicer
        print("i=",i, f[coeff])#print("k=",k,"j=",j,"i=",i, nice_f)
        print("")


##########################################################################

#Functions for the curves y^2=x^(2^m)-1

def gammamatrix(m):
    gammasList=[I] #corresponds to the elliptic curve factor
    for j in [2..m-1]:
        gammasList.append(gammamatrix_x(j))
    return block_diagonal_matrix(gammasList)

def gammamatrixc(m): #the order is the order of gammac for the largest factor = 2^(m-1)
    gammasList=[I] #corresponds to the elliptic curve factor
    for j in [2..m-1]:
        gammasList.append(gammacx(j))
    return block_diagonal_matrix(gammasList)

def gammamatrixcOTHER(m): #the order is the order of gammac for the largest factor = 2^(m-1)
    gammasList=[I] #corresponds to the elliptic curve factor
    for j in [2..m-1]:
        gammasList.append(gammacxOLD(j))
    return block_diagonal_matrix(gammasList)

def Jmatrix(m):
    g = Integer(2^(m-1)-1) #genus of the curve
    blocks=[]
    for i in [1..g]:
        blocks.append(J)
    return block_diagonal_matrix(blocks)


def ID(m):
    F.<z>=CyclotomicField(2^m)
    g=Integer(2^(m-1)-1)
    P=PolynomialRing(F,g,'u')
    IDblocks=[]
    ECblock=[]
    ECblock.append(diagonal_matrix([P.gens()[1],(P.gens()[1])^(-1)]))
    #print(IDblocks)
    count=2
    for d in [2..m-1]:
        IDsmallgenus=[]
        dgenus=Integer(2^(d-1))
        for i in range(Integer(dgenus/2)):
            IDsmallgenus.append(diagonal_matrix([P.gens()[count],(P.gens()[count])^(-1)]))
            count+=1
        IDsmallgenus+=IDsmallgenus
        IDblocks+=IDsmallgenus
    return block_diagonal_matrix(ECblock+IDblocks)


#Checking the order of gamma in the component group
def GammaOrder(m):
    gamma=gammamatrix(m)
    maybeorder=Integer(euler_phi(2^(m))/2)
    g = Integer(2^(m-1)-1) #genus of the curve
    for n in [1..maybeorder]:
        if gamma^n==identity_matrix(2*g) or gamma^n==-identity_matrix(2*g):
            return n
    return 0

#This prints the characteristic polynomials for a given m
def charpoly(m, coeff=1):
    print("Here are the characteristic polynomials for y^2=x^",2^m,"-1")
    print("To save space, it only prints the a_n coefficient, with default a_1")
    print("")
    G0=IDold(m)
    gamma=gammamatrix(m)
    gammaorder=GammaOrder(m) #this is the order of gamma!
    Jm=Jmatrix(m) #the order of this in ST/ST^0 is always 2
    for j in [0..gammaorder-1]:
        for k in [0..1]:
            f=(G0*gamma^j*Jm^k).characteristic_polynomial()
            if f[coeff]!=0:
                print("j=",j,"k=",k, f[coeff])
                print("")


def charpolyWITHC(m,coeff=1):
    print("Here are the characteristic polynomials for y^2=x^",2^m,"-c")
    print("To save space, it only prints the a_n coefficient, with default a_1")
    print("")
    G0=IDold(m)
    #G0=ID(m)
    gamma=gammamatrix(m)
    gammac=gammamatrixc(m)
    gammacorder=2^(m-1)
    gammaorder=GammaOrder(m) #this is the order of gamma!
    Jm=Jmatrix(m) #the order of this in ST/ST^0 is always 2
    for i in [0..gammacorder-1]:
        for j in [0..gammaorder-1]:
            for k in [0..1]:
                f=(G0*gamma^j*Jm^k*gammac^i).characteristic_polynomial()
                #nice_f=(f(SR('x'))).expand().collect(x) #this makes the coefficients look a little nicer
                if f[n]!=0:
                    print("i=",i,"j=",j,"k=",k, f[coeff])
                    #print("i=",i,"j=",j,"k=",k, nice_f)
                    print("")





