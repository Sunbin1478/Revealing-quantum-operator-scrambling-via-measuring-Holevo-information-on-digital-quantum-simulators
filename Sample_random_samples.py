import random
import numpy as np
import projectq
from projectq.ops import C, H, X, Rx, Ry, T, Measure, TimeEvolution, QubitOperator, All, CZ
from scanf import scanf
import time
import json

eng = projectq.MainEngine()

def one_pauli(symbol, i):
    return ''.join([symbol, str(i)])

def two_pauli(s1, s2, i1, i2):
    return ' '.join([one_pauli(s1, i1), one_pauli(s2, i2)])

def hamiltonian_ZZ(J, L, PC):
    ham_ZZ = 0 * QubitOperator(one_pauli('Z', 0))
    for i in range(L-1):
        ham_ZZ += J * QubitOperator(two_pauli('Z', 'Z', i, i + 1))
    if PC == 'pc':
        ham_ZZ += J * QubitOperator(two_pauli('Z', 'Z', L - 1, 0))
    return ham_ZZ

def hamiltonian_Z(hz, L):
    ham_Z = 0 * QubitOperator(one_pauli('Z', 0))
    for i in range(L):
        ham_Z += hz * QubitOperator(one_pauli('Z', i))
    return ham_Z

def hamiltonian_X(hx, L):
    ham_X = 0 * QubitOperator(one_pauli('X', 0))
    for i in range(L):
        ham_X += hx * QubitOperator(one_pauli('X', i))
    return ham_X

def hamiltonian(J, hx, hz, L, PC):
    ham = hamiltonian_ZZ(J, L, PC) + hamiltonian_X(hx, L) + hamiltonian_Z(hz, L)
    return ham

def bellpair(q1, q2):
    H | q1
    C(X, 1) | (q1, q2)
    return q1, q2

def ibellpair(q1, q2):
    C(X, 1) | (q1, q2)
    H | q1
    return q1, q2

def even_config(target_qubits, cycle_number, applied_gates, gate_records):
    actions = [(Rx, np.pi / 2), (Ry, np.pi / 2), (T, None)]
    for i in range(L - 1):
        gate_constructor, angle = random.choice(actions)
        while gate_constructor == applied_gates[i][0]:
            gate_constructor, angle = random.choice(actions)
        if gate_constructor == T:
            T | target_qubits[i]
            gate_records.append((cycle_number, i, 'T', angle))
        else:
            gate = gate_constructor(angle)
            gate | target_qubits[i]
            gate_records.append((cycle_number, i, gate_constructor.__name__, angle))
        applied_gates[i] = (gate_constructor, angle)

    for i in range(0, len(target_qubits) - 1, 2):
        CZ | (target_qubits[i], target_qubits[i + 1])
        gate_records.append((cycle_number, i, 'CZ', None))

    return target_qubits, applied_gates

def odd_config(target_qubits, cycle_number, applied_gates, gate_records):
    actions = [(Rx, np.pi / 2), (Ry, np.pi / 2), (T, None)]
    for i in range(L - 1):
        gate_constructor, angle = random.choice(actions)
        while gate_constructor == applied_gates[i][0]:
            gate_constructor, angle = random.choice(actions)
        if gate_constructor == T:
            T | target_qubits[i]
            gate_records.append((cycle_number, i, 'T', angle))
        else:
            gate = gate_constructor(angle)
            gate | target_qubits[i]
            gate_records.append((cycle_number, i, gate_constructor.__name__, angle))
        applied_gates[i] = (gate_constructor, angle)

    for i in range(1, len(target_qubits) - 1, 2):
        CZ | (target_qubits[i], target_qubits[i + 1])
        gate_records.append((cycle_number, i, 'CZ', None))

    return target_qubits, applied_gates


def transform(target_qubits, num_cycles):
    applied_gates = [(None, None) for _ in range(len(target_qubits))]
    gate_records = []

    for cycle_number in range(0, num_cycles):
        if cycle_number % 2 == 0:
            target_qubits, applied_gates = even_config(target_qubits, cycle_number, applied_gates, gate_records)
        else:
            target_qubits, applied_gates = odd_config(target_qubits, cycle_number, applied_gates, gate_records)

    return target_qubits, gate_records

def aqccpt(J, hx, hz, L, Time, PC, I):
    qubits = eng.allocate_qureg(L + 1)
    qubits[I], qubits[L] = bellpair(qubits[I], qubits[L])
    target_qubits = [qubits[i] for i in range(L) if i != I]
    target_qubits, gate_records = transform(target_qubits, num_cycles)
    TimeEvolution((-Time), hamiltonian(J, hx, hz, L, PC)) | qubits
    X | qubits[(L // 2)]
    TimeEvolution(Time, hamiltonian(J, hx, hz, L, PC)) | qubits

    qubits[I], qubits[L] = ibellpair(qubits[I], qubits[L])
    eng.flush()
    return qubits, gate_records

def add_zero_at_position(binary_string, n):
    new_binary_string = binary_string[:n] + '0' + binary_string[n:]
    return new_binary_string

def add_one_at_position(binary_string, n):
    new_binary_string = binary_string[:n] + '1' + binary_string[n:]
    return new_binary_string

def qubit(L, n, I):
    b = str(bin(n))
    Ln = len(b)
    qb = ''
    for i in range(Ln - 2):
        qb = qb + b[Ln - i - 1]
    for i in range(L - 1 - Ln + 2):
        qb = qb + '0'
    qb = add_zero_at_position(qb, I)
    qb = qb + '0'
    return qb

def oplbar(state, L, I):
    a = 0
    for i in range(pow(2, L - 1)):
        prob = eng.backend.get_probability(qubit(L, i, I), state)
        a += prob
    a = 1 - a
    All(Measure) | qubits
    return a

time_start = time.time()
J = 1
print("请输入演化的总格点数目,演化最大时间,横轴磁场强度,纵轴磁场强度,演化计算间隔,是否成环")
L, Tmax, hx, hz, deltat, Pc = scanf("%d %d %f %f %f %d")
m = 1000
opsize = np.zeros(m)
num_cycles = 8
data1 = np.zeros((int(Tmax / deltat) + 1, 1))

gate_records_all = []

for k in range(int(m)):
    time_end = time.time()
    print('计算演化第', k, '次花费的时间', time_end - time_start)
    data = np.zeros((int(Tmax / deltat) + 1, 1))

    I = L // 2

    Time = 1
    qubits, gate_records = aqccpt(J, hx, hz, L, Time, Pc, I)
    bx = oplbar(qubits, L, I)
    opsize[k] = bx
    gate_records_all.append(gate_records)

# 剔除相同的数据
unique_opsize = []
unique_gate_records = []
for i in range(m):
    is_unique = True
    for j in range(i):
        if abs(opsize[i] - opsize[j]) < 0.00001:
            is_unique = False
            break
    if is_unique:
        unique_opsize.append(opsize[i])
        unique_gate_records.append(gate_records_all[i])
with open(f"gate_records_N={L}_hz={hz}.json", "w") as json_file:
    json.dump(unique_gate_records, json_file)
