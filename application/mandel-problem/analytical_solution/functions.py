import math

###################################################
def nthRoot(n,c):
    """This function returns the solution of transcendental equation in form of
       tanx = cx

    see answer in https://math.stackexchange.com/questions/18718/solution-of-tanx-x

    Args:
        n : the n-th root of equation
        c : the coefficient before x

    Returns:
        the n-th root of equation
    """
    Rn = (n+0.5)*math.pi
    return Rn - 1/c/Rn-(3*c-1)/3/math.pow(c,3)/math.pow(Rn,3) - (30*c*c-20*c+3)/15/math.pow(c,5)/math.pow(Rn,5)

def getBetaValue(A1_2_,n=50):
    """return the list of beta values of

    Args:
        A1_2_ : [A1/A2]
        n (int, optional): [the required n-th]. Defaults to 50.

    Returns:
        [list]: [the value of beta]
    """
    return [nthRoot(i,A1_2_) for i in range(n)]