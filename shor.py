from pyqpanda import *

import random
import math
import numpy as np
from typing import List

from fractions import Fraction


class InitQMachine:
    """ Init Quantum machine and pre-alloc needed qubits/cbits

    Instance used as a context that containes qubits/cbits and related info.

    attributes:
        machine: init by init_quantum_machine()
        qubits: alloced qubits(by qAlloc_many)
        cbits: alloced classical bits(by cAlloc_many)
    
    """
    def __init__(self, qubitsCount, cbitsCount, machineType = QMachineType.CPU):
        self.machine = init_quantum_machine(machineType)
        
        self.qubits = self.machine.qAlloc_many(qubitsCount)
        self.cbits = self.machine.cAlloc_many(cbitsCount)

        self.phase_bits_qubits = cbitsCount
        self.function_value_qubits = qubitsCount - cbitsCount
        self.total_qubits = qubitsCount
        
        print(f'Init Quantum Machine with qubits:[{qubitsCount}] / cbits:[{cbitsCount}] Successfully')
    
    def __del__(self):
        destroy_quantum_machine(self.machine)


class QuantumCircuit:
    """ support Class about QCircuit that futher simplify the quantum circuit's construction.

    Mainly used in build the transform function's circuit(second step during period-finding)
    thus implements the simplified calling about [H,X,CNOT,Toffoli] and BARRIER
    which return self thus achieved chained calls.

    use QuantumCircuit.circuit to get the corresponding QCircuit instance
    
    """
    def __init__(self, qubits, circuit = None, *, expand = 'no'):
        self.qubits = qubits
        self.circuit = circuit if circuit != None else create_empty_circuit()
        
        # QuantumCircuit.expand's value should be one of: ['no', 'expand', 'approximate']
        #
        # which means the strategy to solve the Toffoli gate in the circuit.
        #   - 'no': default. keep the Toffoli gate without expansion.
        #   - 'expand': expand(replace) the Toffoli gate appearing in the circuit to 15-gates expansion form
        #   - 'approximate': will expand the Toffoli gate into a relative phase Toffoli gates
        # 
        # self.expand = expand # default 'no'
        self.expand = 'approximate'
    
    # Single qubit gates
    
    def x(self, index):
        self.circuit << X(self.qubits[index])

        return self
    
    def h(self, index):
        self.circuit << H(self.qubits[index])

        return self
    
    # Multi qubits gates
    
    def cnot(self, control, target):
        if not isinstance(control, list):
            self.circuit << CNOT(self.qubits[control], self.qubits[target])
        else:
            if control[0] not in ['black', 'white']:
                raise ValueError('control pattern is neither "black" or "white"')
            
            if control[0] == 'white':
                # 'white' dot => this control qubit expects an input: 0
                self.x(control[1])
                self.circuit << CNOT(self.qubits[control[1]], self.qubits[target])
                self.x(control[1])
            else:
                # 'black' dot => this control qubit expects an input: 1 (default)
                self.circuit << CNOT(self.qubits[control[1]], self.qubits[target])
        
        return self
    
    def tof(self, control1, control2, target):
        if not isinstance(control1, list):
            control1 = ['black', control1]
        if not isinstance(control2, list):
            control2 = ['black', control2]
        
        if (control1[0] not in ['black', 'white']) or (control2[0] not in ['black', 'white']):
            raise ValueError('control pattern is neither "black" or "white"')
            
        # 'white' dot => this control qubit expects an input: 0
        if control1[0] == 'white':
            self.x(control1[1])
        if control2[0] == 'white':
            self.x(control2[1])
            
        q1 = self.qubits[control1[1]]
        q2 = self.qubits[control2[1]]
        qt = self.qubits[target]
        
        # expand strategies have: ['no'(default), 'expand', 'approximate']
        #
        if self.expand == 'no':
            self.circuit << Toffoli(q1, q2, qt)
        elif self.expand == 'expand':
            # expand Toffoli gate into a traditional implements using 15 gates.
            self.circuit << H(qt) << CNOT(q2, qt) \
                << T(qt).dagger() << CNOT(q1, qt) \
                << T(qt) << CNOT(q2, qt) \
                << T(qt).dagger() << CNOT(q1, qt) \
                << T(qt) << H(qt) \
                << T(q2) << CNOT(q1, q2) \
                << T(q2).dagger() << T(q1) << CNOT(q1, q2)
        elif self.expand == 'approximate':
            # expand the Toffoli gate into a relative phase Toffoli gates
            self.circuit << RY(qt, np.pi / 4) << CNOT(q2, qt) \
                 << RY(qt, np.pi / 4) << CNOT(q1, qt) \
                 << RY(qt, -np.pi / 4) << CNOT(q2, qt) \
                 << RY(qt, -np.pi / 4)
            
        if control1[0] == 'white':
            self.x(control1[1])
        if control2[0] == 'white':
            self.x(control2[1])
        
        return self
    
    # another
    
    def barrier(self, targets = None):
        targets = targets if targets != None else self.qubits
        self.circuit << BARRIER(targets)

        return self


def create_transform_circuit(a, N, ctx):
    """ create a quantum function circuit applied in period-finding's second step 

    In this experiment, we just consider the case that N=21
    thus call with N != 21 will just raise an Error. 

    what's more, {a} should be pre-checked which means {a} can
    just be one of [2,4,5,8,10,11,13,16,17,19,20]

    cites: https://arxiv.org/pdf/1202.6614.pdf

    Args:
        a: a pre-seleted number to run shor-algorithm
        N: A composite number need to factoring
        ctx: Context about global quantum machine(qubits/cbits etc.)
    -------
    Returns:
        constructed transform function's QCircuit

    """
    if N != 21:
        raise ValueError(f'transform circuit for N={N} haven\'t implemented')
    
    qubits = ctx.qubits
    transform_circuit = QuantumCircuit(qubits) # default: expand: 'no'

    # control pattern 
    #   => BLACK: control qubit expected an input 1(default) / WHITE: expected 0
    BLACK = 'black'; WHITE = 'white'

    if a == 2:
        # remember should return QuantumCircuit.circuit => QCircuit
        # the same below
        # 
        return transform_circuit \
            .cnot(2, 1) \
            .tof(0, 1, 6).tof([WHITE, 0], [WHITE, 1], 3) \
            .cnot(2, 1) \
            .tof([WHITE, 1], [BLACK, 2], 5).tof([WHITE, 0], [BLACK, 5], 7) \
            .tof(0, 5, 3).tof([WHITE, 1], [BLACK, 2], 5) \
            .tof([BLACK, 1], [WHITE, 2], 4).tof([WHITE, 0], [BLACK, 4], 5) \
            .cnot(0, 4).cnot(5, 4) \
            .barrier().circuit
    elif a == 4:
        return transform_circuit \
            .cnot(2, 1) \
            .tof([BLACK, 0], [WHITE, 1], 5).tof([WHITE, 0], [WHITE, 1], 3) \
            .cnot(2, 1) \
            .tof([WHITE, 1], [BLACK, 2], 6).tof([WHITE, 0], [BLACK, 6], 5).tof(0, 6, 7) \
            .cnot(2, 1).cnot(1, 6).cnot(2, 1) \
            .tof([WHITE, 0], [BLACK, 6], 7).tof(0, 6, 3).tof([BLACK, 1], [WHITE, 2], 6) \
            .barrier().circuit
    elif a == 5:
        return transform_circuit \
            .tof([WHITE, 1], [BLACK, 2], 7).tof(0, 7, 3) \
            .tof([BLACK, 1], [WHITE, 2], 5).tof(0, 5, 7) \
            .cnot(2, 1) \
            .tof([BLACK, 0], [WHITE, 1], 5) \
            .cnot([WHITE, 1], 3).cnot(2, 1) \
            .barrier().circuit
    elif a == 8:
        return transform_circuit \
            .cnot(0, 6).cnot([WHITE, 0], 3) \
            .barrier().circuit
    elif a == 10:
        return transform_circuit \
            .tof([WHITE, 1], [BLACK, 2], 3).tof([BLACK, 0], [WHITE, 3], 6).tof(0, 3, 7) \
            .cnot(2, 1).cnot(6, 4).cnot(1, 3).cnot(2, 1) \
            .tof([WHITE, 0], [BLACK, 3], 7).tof([BLACK, 1], [WHITE, 2], 3) \
            .cnot(7, 5).cnot(2, 1).cnot(1, 5) \
            .tof(0, 1, 4) \
            .cnot(1, 0).cnot([WHITE, 0], 3).cnot(1, 0).cnot(2, 1) \
            .barrier().circuit
    elif a == 11:
        return transform_circuit \
            .tof([BLACK, 1], [WHITE, 2], 3).tof([WHITE, 0], [BLACK, 3], 7).tof([BLACK, 0], [WHITE, 3], 4) \
            .cnot(2, 1).cnot(1, 3).cnot(2, 1) \
            .tof([WHITE, 0], [BLACK, 3], 5).tof([BLACK, 0], [WHITE, 3], 6).tof([WHITE, 1], [BLACK, 2], 3) \
            .cnot(2, 1).cnot([WHITE, 1], 3).cnot(2, 1) \
            .barrier().circuit
    elif a == 13:
        return transform_circuit \
            .cnot(0, 6).cnot(0, 5).x(3) \
            .barrier().circuit
    elif a == 16:
        return transform_circuit \
            .tof([WHITE, 1], [BLACK, 2], 6) \
            .cnot(2, 1).cnot(6, 4).cnot(1, 4) \
            .tof([WHITE, 0], [BLACK, 6], 7).tof([BLACK, 0], [WHITE, 1], 7) \
            .tof([WHITE, 0], [BLACK, 4], 5).tof(0, 6, 5) \
            .tof([WHITE, 0], [WHITE, 1], 3).tof(0, 4, 3) \
            .cnot(1, 4).cnot(2, 1).cnot(6, 4) \
            .tof([WHITE, 1], [BLACK, 2], 6) \
            .barrier().circuit
    elif a == 17:
        return transform_circuit \
            .tof([BLACK, 1], [WHITE, 2], 6) \
            .cnot(2, 1).cnot(6, 4).cnot(1, 4) \
            .tof([BLACK, 0], [WHITE, 1], 7) \
            .cnot([WHITE, 1], 3).cnot(4, 5).cnot(6, 7) \
            .tof(0, 6, 5).tof(0, 4, 3) \
            .cnot(1, 4).cnot(2, 1).cnot(6, 4) \
            .tof([BLACK, 1], [WHITE, 2], 6) \
            .barrier().circuit
    elif a == 19:
        return transform_circuit \
            .tof([BLACK, 1], [WHITE, 2], 5).tof(0, 5, 3) \
            .cnot(2, 1).cnot(3, 4).cnot(5, 6) \
            .cnot(0, 4).cnot([WHITE, 1], 3).cnot(1, 6) \
            .tof([WHITE, 0], [BLACK, 6], 7).tof([BLACK, 0], [WHITE, 1], 7) \
            .cnot(1, 6) \
            .tof(0, 1, 6) \
            .cnot(2, 1).cnot(5, 6) \
            .barrier().circuit
    elif a == 20:
        return transform_circuit \
            .cnot(0, 7).cnot(0, 5).cnot([WHITE, 0], 3) \
            .barrier().circuit
    else:
        raise ValueError(f'{a} is not a proper option for N=21!')


def create_program(a, N, ctx):
    """ create qpanda QProg (shor-algorithm quantum circuits)

    input a is pre-checked and guaranted to be proper.

    Args:
        a: a pre-seleted number to run shor-algorithm
        N: A composite number need to factoring
        ctx: Context about global quantum machine(qubits/cbits etc.)
    -------
    Returns:
        constructed QProg

    """
    # Step 0. prepare related environment info
    qubits = ctx.qubits
    cbits = ctx.cbits
    phase_bits_qubits = ctx.phase_bits_qubits
    function_value_qubits = ctx.function_value_qubits

    # Step 1. create superposition of states
    #
    # by applying Hadamard gates
    #
    hadamard_circuit = create_empty_circuit()

    for i in range(phase_bits_qubits):
        hadamard_circuit << H(qubits[i])
        
    hadamard_circuit << BARRIER(qubits)

    # Step 2. Implement a unitary transform function
    #
    transform_circuit = create_transform_circuit(a, N, ctx)

    # Step 3. perform a inverse quantum Fourier tranform
    # 
    qft_dagger_circuit = create_empty_circuit()

    for i in range(phase_bits_qubits - 1):
        qft_dagger_circuit << H(qubits[i])
        
        for j in range(i + 1, phase_bits_qubits):
            qft_dagger_circuit << U1(qubits[i], np.pi / (2 ** (j - i))).control(qubits[j])
        
        qft_dagger_circuit << BARRIER(qubits)

    qft_dagger_circuit << H(qubits[phase_bits_qubits - 1])

    # Step 4. build full circuit program
    #
    prog = create_empty_qprog()
    prog << hadamard_circuit << transform_circuit << qft_dagger_circuit

    return prog


def period_finding(a, N, ctx):
    """ period_finding subroutine called in shor-algorithm

    General process:
        init circuit => measure(get phase)
            => calculate period r (by continued fraction expansion)
            => validate period (a^r ≡ 1 (mod N)) thus return.

    Args:
        a: a pre-seleted number to run shor-algorithm
        N: A composite number need to factoring
        ctx: Context about global quantum machine(qubits/cbits etc.)
    -------
    Returns:
        ok: Boolean mark means success(true) or failure
        period: a proper period.
    
    """
    qubits = ctx.qubits
    cbits = ctx.cbits
    phase_bits_qubits = ctx.phase_bits_qubits

    # Step 1. init period-finding QProg corresponding to a & N
    #
    # while actually in this experiment just implement the case N=21
    # 
    prog = create_program(a, N, ctx)

    # Step 2. Measure it (first phase_bits_qubits qubits) => get phase
    #
    # attention the reading order: qubits[i] => cbits[phase_bits_qubits - i - 1]
    # 
    for i in range(phase_bits_qubits):
        prog << Measure(qubits[i], cbits[phase_bits_qubits - i - 1])

    result = run_with_configuration(prog, cbits, 1)
    print(f'  result: {result}') # like {"101": 1}

    # Convert the reading bits to phase
    #
    # eg. {"101": 1} => ["101"] => 5 => 5 / (2^3) = 0.625
    # 
    phase = int(list(result)[0], 2) / (2 ** phase_bits_qubits)
    print(f'   - corresponding phase: {phase}')

    # Step 3. calculate period r (by continued fraction expansion)
    # 
    # eg. a = 4, phi = 0.625 gives the convergents: {0, 1, 1/2, 2/3, 5/8}
    # while r = 8 is invalid (4^8 mod 21 = 16 !== 1)
    # actually in this case only r = 3 is valid (4^3 mod 21 == 1) (https://arxiv.org/pdf/2103.13855.pdf)
    # 
    # Shor algorithm is designed to work for even orders only,
    # However, for certain choices of square coprime x and odd order, the algorithm works.
    # https://arxiv.org/pdf/1310.6446v2.pdf and https://www.nature.com/articles/nphoton.2012.259 point out it. 
    # 
    # a simple validation way is checking r is even or chosen [a] itself 
    # is a perfect square(https://arxiv.org/pdf/2103.13855.pdf)
    # 
    # we'll use Fraction module and using limit_denominator(limit) to gate r
    # according to the above discussion we'll narrow the limit unitl getting a valid r
    # or narrow to limit = 0 which means fail to find period. 
    #
    limit = N

    while True:
        frac = Fraction(phase).limit_denominator(limit)
        r = frac.denominator

        # simply check period
        if (a ** r) % N == 1:
            break
        else:
            # narrow limit to calculate new period
            limit = r - 1
            if limit <= 0:
                print(f'\n  Rewrite to fraction: {frac}, find period failed')

                return False, None

    print(f'\n  Rewrite to fraction: {frac}, thus r = {r}')

    return True, r


def calculate_factors(a, N, r, ctx):
    """ calculate factors based on calculated period r.

    Possible case: r is even => gcd(a ** (r // 2) - 1, N) ... 
        or r is odd but {a} itself is perfect square => gcd(a ** (r / 2) - 1, N)

    Args:
        a: a pre-seleted number to run shor-algorithm
        N: A composite number need to factoring
        r: calculated valid period r
        ctx: Context about global quantum machine(qubits/cbits etc.)
    -------
    Returns:
        ok: Boolean mark means success(true) or failure
        factors: a pair of factors(should be sorted) when success or empty [] when failed.
     
    """
    # According to Shor algorithm, calculate gcd(a^(r/2) - 1, N) and gcd(a^(r/2) + 1, N)
    # 
    if int(math.sqrt(a)) ** 2 == a:
        # a itself is a perfect square thus Shor still works
        # cite: https://arxiv.org/pdf/2103.13855.pdf
        # 
        guesses = [math.gcd(int(a ** (r / 2)) - 1, N), math.gcd(int(a ** (r / 2)) + 1, N)]
    else:
        guesses = [math.gcd(a ** (r // 2) - 1, N), math.gcd(a ** (r // 2) + 1, N)]
    
    print(f'  calculate final guesses: {guesses}')
    
    # need to check the calculated guesses numbers.
    # 
    factors = []
    for num in guesses:
        if num not in [1, N]:
            factors.append(num)
            print(f'[Find non-trivial factor: {num}]')
    
    if len(factors) == 0:
        print('[Failed]')

        return False, []
    elif len(factors) == 1:
        # may just find one factor and 1
        # calculate another
        # 
        factors.append(N // factors[0])
    
    return True, sorted(factors)


def shor_alg(N, *, a = None, ctx):
    """ Shor algorithm using to factoring a composite number N

    Args:
        N: A composite number need to factoring
        a: a pre-seleted number to run shor-algorithm
        ctx: Context about global quantum machine(qubits/cbits etc.)
    -------
    Returns:
        A pair of non-trivial factors of N (factors should be sorted)
        return {None} if failed.

    """

    # Step 0. If N is even, 2 is the factor
    if N % 2 == 0:
        return [2, N // 2]
    
    # Step 1. Randomly choose a number 1 < a < N
    if a == None:
        a = random.randint(2, N - 1)
        print(f'Randomly choose a = {a}\n')

    # Step 2. Check the randomly choosed number a
    # 
    # compute K = gcd(a, N), 
    # if K != 1, then K is thus a non-trivial factor of N
    # algorithm finished with returned [K, N / K]
    # 
    K = math.gcd(a, N)
    if K != 1:
        # thus K is one of a factor of N
        print(f' - gcd({a}, {N}) = {K}! {K} is a factor of {N}')
        print('\nShor Algorithm Finished!')

        return sorted([K, N // K])

    # Step 3. call quantum period-finding subroutine to find period r
    #
    # should check r that a^r ≡ 1 (mod N) (checked by subroutine function)
    # after getting a proper period r, calculate factors and return.
    #
    # because each time the running result will affected by nondeterministic measurements
    # (during period-finding, will measure the phase result to calculate period) 
    # for selected {a}, will try {MAX_ATTEMPT} to calculate period repeatly
    # 
    MAX_ATTEMPT = 20
    attempt = 1
    
    while True:
        print(f'Execute shor algorithm - attempt time: [{attempt}]')

        # call period-finding subroutine
        valid, r = period_finding(a, N, ctx)
        if valid:
            # valid period, and then calculate factors:
            ok, factors = calculate_factors(a, N, r, ctx)
            if ok:
                print(f'\nFactors of {N} are: {factors} (total run times: {attempt})')
                print('\nShor Algorithm Finished!')
                
                # reutrned factors should be sorted(by subroutine)
                return factors
            
        attempt += 1
        if attempt > MAX_ATTEMPT:
            print('\nShor Algorithm Finished [FAIL]')
            
            return None
        
        print('\n' + '-' * 36 + '\n')
    

def solution() -> List[int]:
    # Step 0. Init qMachine and qubits/cbits
    phase_bits_qubits = 3
    function_value_qubits = 5

    total_qubits = phase_bits_qubits + function_value_qubits

    ## init quantum machine
    ## related info(qubits/cbits and counts etc.) packed into a context object
    # 
    ctx = InitQMachine(total_qubits, phase_bits_qubits)

    # Cause the shor algorithm might fail => attempt some rounds
    MAX_ROUND = 8
    
    N = 21

    for round in range(MAX_ROUND):
        print(f'Attempt call shor-algorithm: round {round + 1}\n')

        res = shor_alg(N, ctx = ctx)

        if res != None:
            return res

    return