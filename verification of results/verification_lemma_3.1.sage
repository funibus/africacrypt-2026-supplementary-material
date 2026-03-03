#############################################
## This sage script verifies the statement ##
## from Lemma 3.1  and Corollary 3.2 in    ##
## the article "On the Properties of       ##
##  of HighBits and LowBits Functions      ##
##  and their Applications"                ##
#############################################

## How to use?
#### Simply run this script in Sage. It should print "All test passed!" if the lemma is correct. If some test did not pass, it will output a pair (q,alpha) for which the lemma is incorrect.

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
    
def Lb(q,alpha,t):
    ## Compute Hb_{q, alpha}(t)
    return P2R(q,alpha,t)[0]
    
def decompose_q(q,alpha):
    ## compute q0 and q1 s.t. q = q0+alpha*q1 with q0 in (-alpha/2,alpha/2]
    q0 = q%alpha
    if q0 > alpha/2:
        q0 = q0-alpha ## q0 in (-alpha/2,alpha/2]
    q1 = (q-q0)/alpha
    return(q0,q1)
       
    
## list of expected values

def allowed_values(q,alpha):
    ## returns the lists of all 6 values that are expected for t^h and t^l
    (q0,q1) = decompose_q(q,alpha)
    res_h = [-1,0,1, -q1-1, -q1, -q1+1]
    res_l = [-alpha,0,alpha, -q0-alpha, -q0, -q0+alpha]
    return (res_h, res_l)


## Test function

def enum_all(q,alpha):
    ## computes t^h and t^l for all pairs (t,r) in [0,q)
    res_h = []
    res_l = []
    for t in range(q):
        for r in range(t,q):
            th = Hb(q, alpha, t+r)-Hb(q, alpha, t)-Hb(q, alpha, r)
            tl = Lb(q, alpha, t+r)-Lb(q, alpha, t)-Lb(q, alpha, r)
            if not th in res_h:
                res_h += [th]
            if not tl in res_l:
                res_l += [tl]
    return (res_h, res_l)

def test_list(list_q_alpha, verbose = False):
    ## for all (q,alpha) in list_q_alpha, compute the lists of all values expected for t^h and t^l, and the lists of all values actually achieved, and check that the second ones are contained in the first ones
    failed = False
    for (q,alpha) in list_q_alpha:
        if verbose:
            print("\n(q,alpha) = ", (q, alpha))
        (expected_th, expected_tl) = allowed_values(q,alpha)
        if verbose:
            print("expected values for th =", expected_th)
            print("expected values for tl =", expected_tl)
        (actual_th, actual_tl) = enum_all(q,alpha)
        if verbose:
            print("actual values observed for th =", actual_th)
            print("actual values observed for tl =", actual_tl)
        for th in actual_th:
            if not th in expected_th:
                print("\n Warning, one t^h was not expected for (q,alpha) = ", (q,alpha))
                failed = True
        for tl in actual_tl:
            if not tl in expected_tl:
                print("\n Warning, one t^l was not expected for (q,alpha) = ", (q,alpha))
                failed = True
    if failed == True:
        print("\nSome test did not pass")
    else:
        print("All tests passed!")
        
### Running the tests

#set_random_seed(42)
list_q_alpha = []
for _ in range(10): ## random pairs of alpha and q
    alpha = ZZ.random_element(5,100)
    q = ZZ.random_element(200,500)
    list_q_alpha += [(q,alpha)]
    
for d in range(3,7): ## random pairs of alpha and q of the form of dilithium
    alpha = 2^d ## alpha is either 8, 16, 32, or 64
    q1 = 3*alpha+1
    q2 = 5*alpha+1
    list_q_alpha += [(q1,alpha), (q2,alpha)]
    
print("Verification of Lemma 3.1")
print("testing 18 random choices of q and alpha (can take a few seconds)...")
test_list(list_q_alpha)

### Random verification of Corollary 3.2

def random_test(q,alpha, nb_tests = 100000):
    failed = False
    (expected_th, expected_tl) = allowed_values(q,alpha)
    for _ in range(nb_tests):
        t = ZZ.random_element(q)
        r = ZZ.random_element(q)
        th = Hb(q, alpha, t+r)-Hb(q, alpha, t)-Hb(q, alpha, r)
        tl = Lb(q, alpha, t+r)-Lb(q, alpha, t)-Lb(q, alpha, r)
        if not th in expected_th:
            print("\n Warning, one t^h was not expected for (t,r) = ", (t,r))
            failed = True
        if not tl in expected_tl:
            print("\n Warning, one t^l was not expected for (t,r) = ", (t,r))
            failed = True
    if failed == True:
        print("\nSome test did not pass")
    else:
        print("All tests passed!")
        
(q,alpha) = (8380417, 2^13)
(expected_th, expected_tl) = allowed_values(q,alpha)
print("\n\nVerification of Corollary 3.2")
print("For (q,alpha) = ", (q, alpha))
print("The expected t^h are: ", expected_th)
print("The expected t^l are: ", expected_tl)
print("Checking on 100000 random values of t and r...")
random_test(q,alpha)
