import random
import numpy as np
import projectq
from projectq.ops import C, H, X, Rx, Ry,T, Measure, TimeEvolution, QubitOperator,  All, CZ
from scanf import scanf
import time
eng = projectq.MainEngine()
def one_pauli(symbol, i):
    return ''.join([symbol, str(i)])

def two_pauli(s1, s2, i1, i2):
    return ' '.join([one_pauli(s1, i1), one_pauli(s2, i2)])

# 定义泡利算符
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
def even_config(target_qubits, cycle_number,applied_gates):
    actions = [(Rx, np.pi / 2), (Ry, np.pi / 2), (T, None)]
    for i in range(L - 1):
        gate_constructor, angle = random.choice(actions)
        while gate_constructor == applied_gates[i][0]:#不与上次循环重复
            gate_constructor, angle = random.choice(actions)
        if gate_constructor == T:
            T | target_qubits[i]
        else:
            gate = gate_constructor(angle)
            gate | target_qubits[i]
        applied_gates[i] = (gate_constructor, angle)
    for i in range(0, len(target_qubits) - 1, 2):
        CZ | (target_qubits[i], target_qubits[i + 1])

    return target_qubits, applied_gates


def odd_config(target_qubits, cycle_number,applied_gates):
    actions = [(Rx, np.pi / 2), (Ry, np.pi / 2), (T,None)]
    for i in range(L-1):
        gate_constructor, angle = random.choice(actions)
        while gate_constructor == applied_gates[i][0]:
            gate_constructor, angle = random.choice(actions)
        if gate_constructor == T:
            T | target_qubits[i]
        else:
            gate = gate_constructor(angle)
            gate | target_qubits[i]
        applied_gates[i] = (gate_constructor, angle)
    for i in range(1, len(target_qubits) - 1, 2):
        CZ | (target_qubits[i], target_qubits[i + 1])

    return target_qubits, applied_gates


def transform(target_qubits, num_cycles):
    applied_gates = [(None, None) for _ in range(len(target_qubits))]

    for cycle_nmuber in range(0, num_cycles):
        if cycle_nmuber % 2 == 0:
            target_qubits, applied_gates  = even_config(target_qubits, cycle_nmuber,applied_gates)
        else:
            target_qubits, applied_gates = odd_config(target_qubits, cycle_nmuber,applied_gates)

    return target_qubits

def bellpair(q1, q2):
    H | q1
    C(X, 1) | (q1, q2)
    return q1, q2

def ibellpair(q1, q2):
    C(X, 1) | (q1, q2)
    H | q1
    return q1, q2

def aqccpt(J, hx, hz, L, Time, PC, I):
    qubits = eng.allocate_qureg(L+1)
    qubits[I], qubits[L] = bellpair(qubits[I], qubits[L])
    target_qubits = [qubits[i] for i in range(L) if i != I]
    target_qubits=transform(target_qubits,num_cycles)
    TimeEvolution((-Time), hamiltonian(J, hx, hz, L, PC)) | qubits
    X | qubits[(L // 2)]
    TimeEvolution(Time, hamiltonian(J, hx, hz, L, PC)) | qubits

    qubits[I], qubits[L] = ibellpair(qubits[I], qubits[L])
    eng.flush()
    return qubits

def add_zero_at_position(binary_string, n)
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
    for i in range(pow(2, L - 1)):#得到结为贝尔态的格点投影全为零的概率
        prob =eng.backend.get_probability(qubit(L, i, I), state)
        a += prob
    return a
def qubitX(L, n, I):#得到总长度为2*L，最后第I位为0和L+I位为1的字符串
    b = str(bin(n))
    Ln = len(b)
    qb = ''
    for i in range(Ln - 2):
        qb = qb + b[Ln - i - 1]
    for i in range(L - 1 - Ln + 2):
        qb = qb + '0'
    qb = add_zero_at_position(qb, I)
    qb = qb + '1'
    return qb


def oplbarX(state, L, I):
    a = 0
    for i in range(pow(2, L - 1)):
        prob =eng.backend.get_probability(qubitX(L, i, I), state)
        a += prob
    return a
def qubitY(L, n, I):
    b = str(bin(n))
    Ln = len(b)
    qb = ''
    for i in range(Ln - 2):
        qb = qb + b[Ln - i - 1]
    for i in range(L - 1 - Ln + 2):
        qb = qb + '0'
    qb = add_one_at_position(qb, I)
    qb = qb + '1'
    return qb

def oplbarY(state, L, I):
    a = 0
    for i in range(pow(2, L - 1)):
        prob =eng.backend.get_probability(qubitY(L, i, I), state)
        a += prob

    return a
def qubitZ(L, n, I):
    b = str(bin(n))
    Ln = len(b)
    qb = ''
    for i in range(Ln - 2):
        qb = qb + b[Ln - i - 1]
    for i in range(L - 1 - Ln + 2):
        qb = qb + '0'
    qb = add_one_at_position(qb, I)
    qb = qb + '0'
    return qb

def oplbarZ(state, L, I):
    a = 0
    for i in range(pow(2, L - 1)):#得到结为贝尔态的格点投影全为零的概率
        prob =eng.backend.get_probability(qubitZ(L, i, I), state)
        a += prob
    return a
time_start = time.time()  # 计时开始
def save_data_to_txt(data, filename='data.txt'):

    with open(filename, 'w') as file:
        for row in data:
            row_str = ', '.join(map(str, row))
            file.write(f'[{row_str}]\n')


J = 1
num_cycles=8
m=100
print("请输入演化的总格点数目,演化最大时间,横轴磁场强度,纵轴磁场强度,演化计算间隔,是否成环")
L, Tmax, hx, hz, deltat, Pc = scanf("%d %d %f %f %f %d")

# 生成一个零矩阵,便于存入数据
time_end = time.time()  # 计时结束

data = np.zeros((int(Tmax / deltat) + 1, int(L//2)+1))
dataX = np.zeros((int(Tmax / deltat) + 1, int(L//2)+1))
dataY = np.zeros((int(Tmax / deltat) + 1, int(L//2)+1))
dataZ = np.zeros((int(Tmax / deltat) + 1, int(L//2)+1))
for j in range(int(Tmax / deltat) + 1):
    Time = j * deltat
    time_end = time.time()  # 计时结束
    print('计算演化第',Time, '时刻花费的时间', time_end - time_start)  # 输出时间
    for k in range(int(m)):
        eng.flush()
        for i in range((L//2)+1):
            I = i
            qubits = aqccpt(J, hx, hz, L, Time, Pc,I)
            bx = oplbar(qubits, L, I)
            Xx = oplbarX(qubits, L, I)
            Yx = oplbarY(qubits, L, I)
            Zx = oplbarZ(qubits, L, I)
            data[j, i] += (bx / m)
            dataX[j, i] += (Xx / m)
            dataY[j, i] += (Yx / m)
            dataZ[j, i] += (Zx / m)
            All(Measure) | qubits  # 便于flush清除比特
            eng.flush()
HolevoI=str('HolevoI')+str(',T=')+str(Tmax)+str(',N=')+str(L)+str(',hz=')+str(hz)+str(',hx=')+str(hx)+str(',σ=')+str("X")
save_data_to_txt(data,filename=HolevoI)
HolevoX=str('HolevoX')+str(',T=')+str(Tmax)+str(',N=')+str(L)+str(',hz=')+str(hz)+str(',hx=')+str(hx)+str(',σ=')+str("X")
save_data_to_txt(dataX,filename=HolevoX)
HolevoY=str('HolevoY')+str(',T=')+str(Tmax)+str(',N=')+str(L)+str(',hz=')+str(hz)+str(',hx=')+str(hx)+str(',σ=')+str("X")
save_data_to_txt(dataY,filename=HolevoY)
HolevoZ=str('HolevoZ')+str(',T=')+str(Tmax)+str(',N=')+str(L)+str(',hz=')+str(hz)+str(',hx=')+str(hx)+str(',σ=')+str("X")
save_data_to_txt(dataZ,filename=HolevoZ)


time_end = time.time()  # 计时结束


