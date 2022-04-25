from config import *

# return only up/down spin part of state vector
def up(x): return x[0:NSITES]
def down(x): return x[NSITES:2*NSITES]

# returns total occupation of a site
def n_site(x, site):
    if up(x)[site] == 1 and down(x)[site] == 1: # double occupied
        return 2
    elif up(x)[site] == 1 or down(x)[site] == 1: # single occupied
        return 1
    return 0 # no occupation

# returns list of binary digits of an integer i, padded on the left with zeros
def toBinary(i, length):
    bin_i = bin(i)[2:]
    padding = length - len(bin_i)
    #print([i, bin_i, padding, np.array(padding*[0] + [int(x) for x in bin_i])])
    return np.array(padding*[0] + [int(x) for x in bin_i])

# generate all possible states of the system and save states to a list
def generate_states(N_electrons):
    states = []
    for i in range(0, 4**NSITES):
        state = toBinary(i, 2*NSITES)
        if sum(state) == N_electrons:
            states += [state]
            #print(state)
    return states

# generate pairs of spin sites with allowed hoppings
hoppings = np.array([[NSITES-1,0] if i==NSITES-1 else [2*NSITES-1,NSITES] if i==2*NSITES-1 else [i,i+1] for i in range(2*NSITES)])
hoppings = np.append(hoppings, np.flip(hoppings, 1), 0)
hoppings = hoppings.tolist()

# hop check between two states
def hop_check(state1, state2):
    state_diff = state1 - state2
    
    if sum(abs(state_diff)) == 2:
        test_hop = np.array([], dtype=int)
        for i, spinsite in enumerate(state_diff):
            if abs(spinsite) == 1:
                test_hop = np.append(test_hop, [i*spinsite])
            test_hop = np.sort(test_hop)
        return abs(test_hop).tolist()
    return [-1, -1]

# get the sign for a possible hop between two states
def hop_sign(test_hop, state1, state2):
    return (-1)**sum(state1[0:test_hop[0]]) * (-1)**sum(state2[0:test_hop[1]])

# generate the Hamiltonian
def generate_Hamiltonian(states, Ui=U):
    dim = len(states)
    H_coulomb, H_hop = np.zeros((dim, dim)), np.zeros((dim, dim))

    for i, state in enumerate(states):
        for j in range(0, 3):            
            if n_site(state, j) == 2:
                H_coulomb[i][i] += Ui
            
    for i, state1 in enumerate(states):
        for j, state2 in enumerate(states):
            test_hop = hop_check(state1, state2)
            
            if test_hop in hoppings:
                H_hop[i, j] += hop_sign(test_hop, state1, state2) * t
    
    return H_coulomb + H_hop

# generate the Sz matrix as described in the task
def calc_Sz(states):
    i = 0
    j = 0
    Sz_a, Sz_b, Sz_c = np.zeros([len(states), len(states)]), np.zeros([len(states), len(states)]), np.zeros([len(states), len(states)])

    while j < NSITES:    
        while i < len(states):
            res = up(states[i])[j] - down(states[i])[j]
            if j == 0:
                Sz_a[i, i] = res
            if j == 1:
                Sz_b[i, i] = res
            if j == 2:
                Sz_c[i, i] = res
            i+=1
        i=0
        j+=1
    return (Sz_a*Sz_b + Sz_b*Sz_c + Sz_a*Sz_c)/3