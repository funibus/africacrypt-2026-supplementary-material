#############################################
## This sage script verifies the statement ##
## from Theorem 3.6 in the article "On     ##
## the Properties of HighBits and LowBits  ##
## Functions and their Applications"       ##
#############################################

## How to use?
#### Simply run this script in Sage. It should print "All test passed!" if the theorem is correct. If some test did not pass, it will output a pair (q,alpha) for which the theorem is incorrect.
#### The number of tests is 40 in total: 10 for each possible interval choice of q0.

## Compute the Hb function from the article

def P2R(q,alpha,t):
    ## Compute Power2Round{q,alpha}(t)
    tq = t%q ## tq in [0,q)
    t0 = tq%alpha
    if t0 > alpha/2:
        t0 = t0-alpha ## t0 in (-alpha/2, alpha/2]
    t1 = (tq-t0)/alpha
    return(t0,t1)
    
def Hb(q,alpha,t):
    ## Compute Hb_{q, alpha}(t)
    return P2R(q,alpha,t)[1]
    
def decompose_q(q,alpha):
    ## compute q0 and q1 s.t. q = q0+alpha*q1 with q0 in (-alpha/2,alpha/2]
    q0 = q%alpha
    if q0 > alpha/2:
        q0 = q0-alpha ## q0 in (-alpha/2,alpha/2]
    q1 = (q-q0)/alpha
    return(q0,q1)
    
def HbS(q,alpha,t,r):
    ## Compute Hb(t+r)
    return Hb(q,alpha,t+r)
    
def HbX(q,alpha,t,r):
    ## Compute Hb((Hb(t)+Hb(r))*alpha)
    t1 = Hb(q,alpha,t)
    r1 = Hb(q,alpha,r)
    return Hb(q, alpha, alpha*(t1+r1))
    
## compute a list of alpha and q's

def random_pair(interval):
    ## interval can be:
    ##    1  -> for q0 in (-alpha/2,-alpha/4)
    ##    2  -> for q0 in [-alpha/4,0]
    ##    3  -> for q0 in (0,alpha/2)
    ##    4  -> for q0 = alpha/2
    ## the algorithm generates a random pair (q,alpha) such that q0 lives in the desired interval
    alpha = ZZ.random_element(5,100)
    if interval == 4 and alpha%2 != 0: ## for interval 4, we need alpha to be even
        alpha = alpha+1
    q = ZZ.random_element(3,10)*alpha
    
    lb1 = floor(-alpha/2+1) ## lower bound of first interval
    lb2 = ceil(-alpha/4)
    ub3 = ceil(alpha/2-1)
    if interval == 1:
         q = q+ZZ.random_element(lb1, lb2)
    if interval == 2:
        q = q+ZZ.random_element(lb2, 1)
    if interval == 3:
        q = q+ZZ.random_element(1,ub3)
    if interval == 4:
        q = q+alpha/2
        
    ## verification
    (q0,q1) = decompose_q(q,alpha)
    if interval == 1:
        assert(-alpha/2 < q0 and q0 < -alpha/4)
    if interval == 2:
        assert(-alpha/4 <= q0 and q0 <= 0)
    if interval == 3:
        assert(0 < q0 and q0 < alpha/2)
    if interval == 4:
        assert(q0 == alpha/2)
    return(q,alpha)
    
## list of expected values

def allowed_values(q,alpha):
    ## returns the list of all 5,6 or 7 values that are expected for S1-X1, depending on the choice of q and alpha
    (q0,q1) = decompose_q(q,alpha)
    res = [-1,0,1]
    if -alpha/2 < q0 and q0 < -alpha/4:
        res += [-q1+1,q1-2,q1-1,q1]
    elif -alpha/4 <= q0 and q0 <= 0:
        res += [-q1+1,q1-1,q1]
    elif 0 < q0 and q0 < alpha/2:
        res += [-q1,-q1+1,q1-1]
    elif q0 == alpha/2:
        res += [-q1,q1]
    return res


## Test function

def enum_all(q,alpha):
    ## computes eps = S1-X1 all pairs (t,r) in [0,q)
    res = {}
    for t in range(q):
        for r in range(t,q):
            eps = HbS(q,alpha,t,r) - HbX(q,alpha,t,r)
            if not eps in res:
                res[eps] = 0
            res[eps] += 1
    return res

def test_list(list_q_alpha, verbose = False):
    ## for all (q,alpha) in list_q_alpha, compute the list of all values expected for eps, and the list of all values actually achieved, and check that the second one is contained in the first one
    failed = False
    for (q,alpha) in list_q_alpha:
        if verbose:
            print("\n(q,alpha) = ", (q, alpha))
        expected_eps = allowed_values(q,alpha)
        if verbose:
            print("expected values for eps =", expected_eps)
        actual_eps = enum_all(q,alpha)
        if verbose:
            print("actual values observed (and their frequency) =", actual_eps)
        for eps in actual_eps:
            if not eps in expected_eps:
                print("\n Warning, one eps was not expected for (q,alpha) = ", (q,alpha))
                failed = True
    if failed == True:
        print("\nSome test did not pass")
    else:
        print("All tests passed!")
        
### Running the tests

#set_random_seed(42)
list_q_alpha = []
for _ in range(10):
    list_q_alpha += [random_pair(1), random_pair(2), random_pair(3), random_pair(4)]
    
print("testing 10 random choices of q and alpha for each possible interval of q0 (can take a few minutes)...")
test_list(list_q_alpha)
    
