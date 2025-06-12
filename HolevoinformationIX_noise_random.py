import random
from qibo import hamiltonians, models
from qibo.symbols import X, Y, Z
import os
import numpy as np
from qibo.models import Circuit
from qibo import gates
import time
import multiprocessing as mp
from functools import partial
os.environ["QIBO_BACKEND"] = "numpy"  # 强制使用 NumPy 后端
def list(p):
    return [("X", p / 3), ("Y", p / 3), ("Z", p / 3)]
NL=7
def bellst(c, n, m, p):
    c.add(gates.H(n))
    c.add(gates.PauliNoiseChannel(n, list(p)))
    c.add(gates.CNOT(n, m))
    c.add(gates.PauliNoiseChannel(n, list(p)))
    c.add(gates.PauliNoiseChannel(m, list(p)))
    return c

def ibellst(c, n, m, p):
    c.add(gates.CNOT(n, m))
    c.add(gates.PauliNoiseChannel(n, list(p)))
    c.add(gates.PauliNoiseChannel(m, list(p)))
    c.add(gates.H(n))
    c.add(gates.PauliNoiseChannel(n, list(p)))
    return c

def RZZ(J, dt, c, q1, q2, p):
    c.add(gates.CNOT(q1, q2))
    c.add(gates.PauliNoiseChannel(q1, list(p)))
    c.add(gates.PauliNoiseChannel(q2, list(p)))
    c.add(gates.RZ(q2, theta=-2 * J * dt))
    c.add(gates.PauliNoiseChannel(q2, list(p)))
    c.add(gates.CNOT(q1, q2))
    c.add(gates.PauliNoiseChannel(q1, list(p)))
    c.add(gates.PauliNoiseChannel(q2, list(p)))
    return c

def Hmix_dt_evl(c, NL, J, hz, hx, dt, p, direction):
    for i in range(int(NL / 2)):
        c = RZZ(J, dt, c, 2 * i, 2 * i + 1, p)
    for i in range(int((NL - 1) / 2)):
        c = RZZ(J, dt, c, 2 * i + 1, 2 * i + 2, p)
    for i in range(NL):
        theta = -2 * hz * dt if direction == 'up' else 2 * hz * dt
        c.add(gates.RZ(i, theta=theta))
        c.add(gates.PauliNoiseChannel(i, list(p)))
    for i in range(NL):
        theta = -2 * hx * dt if direction == 'up' else 2 * hx * dt
        c.add(gates.RX(i, theta=theta))
        c.add(gates.PauliNoiseChannel(i, list(p)))
    return c
def even_config(circuit, target_qubits,  applied_gates,p):
    actions = [(gates.RX, np.pi / 2), (gates.RY, np.pi / 2), (gates.T, None)]
    for i in range(NL-1):
        gate_constructor, angle = random.choice(actions)
        while gate_constructor == applied_gates[i][0]:
            gate_constructor, angle = random.choice(actions)
        if gate_constructor == gates.T:
            circuit.add(gates.T(target_qubits[i]))
            circuit.add(gates.PauliNoiseChannel(target_qubits[i], list(p)))

        else:
            gate = gate_constructor(target_qubits[i], angle)
            circuit.add(gate)
            circuit.add(gates.PauliNoiseChannel(target_qubits[i], list(p)))
        applied_gates[i] = (gate_constructor, angle)
    for i in range(0, len(target_qubits) - 1, 2):
        circuit.add(gates.CZ(target_qubits[i], target_qubits[i + 1]))
        circuit.add(gates.PauliNoiseChannel(target_qubits[i], list(p)))
        circuit.add(gates.PauliNoiseChannel(target_qubits[i + 1], list(p)))
    return applied_gates

def odd_config(circuit, target_qubits,  applied_gates,p):
    actions = [(gates.RX, np.pi / 2), (gates.RY, np.pi / 2), (gates.T, None)]
    for i in range(NL-1):
        gate_constructor, angle = random.choice(actions)
        while gate_constructor == applied_gates[i][0]:
            gate_constructor, angle = random.choice(actions)
        if gate_constructor == gates.T:
            circuit.add(gates.T(target_qubits[i]))
            circuit.add(gates.PauliNoiseChannel(target_qubits[i], list(p)))
        else:
            gate = gate_constructor(target_qubits[i], angle)
            circuit.add(gate)
            circuit.add(gates.PauliNoiseChannel(target_qubits[i], list(p)))

        applied_gates[i] = (gate_constructor, angle)
    for i in range(1, len(target_qubits) - 1, 2):
        circuit.add(gates.CZ(target_qubits[i], target_qubits[i + 1]))
        circuit.add(gates.PauliNoiseChannel(target_qubits[i], list(p)))
        circuit.add(gates.PauliNoiseChannel(target_qubits[i+1], list(p)))

    return applied_gates

def transform(circuit,target_qubits, num_cycles, p):
    applied_gates = [(None, None) for _ in range(len(target_qubits))]
    gate_records = []
    for cycle_number in range(num_cycles):
        if cycle_number % 2 == 0:
            applied_gates = even_config(circuit, target_qubits,  applied_gates,p)
        else:
            applied_gates = odd_config(circuit, target_qubits,  applied_gates,p)
    return  circuit
num_cycles=8
def time_evl_get_Lbar(NL, J, hx, hz, T, ci, p, I, gate_type):
    step = max(1, int(T / 0.1))
    dt = T / step
    p = ci * p
    c = models.Circuit(NL + 1, density_matrix=True)
    c = bellst(c, I, NL, p)
    target_qubits = [j for j in range(NL) if j != I]
    c = transform(c, target_qubits, num_cycles, p)
    for _ in range(step):
        c = Hmix_dt_evl(c, NL, J, hz, hx, dt, p, direction='down')

    if gate_type == 'X':
        c.add(gates.X(int(0)))
    else:
        c.add(gates.Y(int(0)))
    c.add(gates.PauliNoiseChannel(int(0), list(p)))

    for _ in range(step):
        c = Hmix_dt_evl(c, NL, J, hz, hx, -dt, p, direction='up')

    c = ibellst(c, I, NL, p)
    result = c()
    return oplbar(result, I, NL)

def entropy(p):
    return -np.sum(p * np.log2(p + 1e-14))  # 添加小量以避免log(0)
def holevo_information(L):
    return 0.5 * (1 + np.log2(2 / (2 - L) * ((1 - L) / (2 - L)) ** (1 - L)))
def oplbar(result, NL, I):
    A = result.state()
    qubit_indices = [I, NL]
    prob_00 = probability_00(A, qubit_indices)
    return prob_00
def probability(state, qubit_indices, target_state):
    num_qubits = int(np.log2(len(state)))
    total_prob = 1.0
    for i in range(len(state)):
        binary_state = f"{i:0{num_qubits}b}"
        if all(binary_state[qubit_indices[j]] == target_state[j] for j in range(len(qubit_indices))):
            total_prob -= np.real(state[i, i])
    return total_prob
def probability_00(state, qubit_indices):
    return probability(state, qubit_indices, '00')


# 新增：计算单个时间步的函数
def compute_time_step(i, NL, J, hx, hz, deltaT, p):
    T = deltaT * (i + 1)
    I = 0
    E0X_00 = 0
    E1X_00 = 0
    E2X_00 = 0
    E3X_00 = 0
    E4X_00 =0
    for w in range(100):
        # 使用相同的量子电路实例进行多次计算
        E0X_00 += time_evl_get_Lbar(NL, J, hx, hz, T, 1, 1e-10, I, 'X')
        E1X_00 += time_evl_get_Lbar(NL, J, hx, hz, T, 1, p, I, 'X')
        E2X_00 += time_evl_get_Lbar(NL, J, hx, hz, T, 2, p, I, 'X')
        E3X_00 += time_evl_get_Lbar(NL, J, hx, hz, T, 3, p, I, 'X')
        E4X_00 += time_evl_get_Lbar(NL, J, hx, hz, T, 4, p, I, 'X')
    E0 =E0X_00/100
    E1 =E1X_00/100
    E2 = (2 * E1X_00 - E2X_00)/100
    E3 = (3 * E1X_00 - 3 * E2X_00 + E3X_00)/100
    E4 = (4 * E1X_00 - 6 * E2X_00 + 4 * E3X_00 - E4X_00)/100

    print(f'完成计算演化第 {T} 时刻')
    return i, E0, E1, E2, E3, E4


# 新增：批量写入文件函数
def write_results(base_dir, deltaT, results, Tmax, NL, hz, p):
    # 准备所有数据
    E0, E1, E2, E3, E4 = zip(*[(r[1], r[2], r[3], r[4], r[5]) for r in sorted(results)])

    # 准备文件数据
    files_data = [
        (f'HtestIX,T={Tmax},N={NL},hz={hz},p={1e-10},noiseless.dat', E0),
        (f'HtestIX,T={Tmax},N={NL},hz={hz},p={p},mtg-n=0.dat', E1),
        (f'HtestIX,T={Tmax},N={NL},hz={hz},p={p},mtg-n=1.dat', E2),
        (f'HtestIX,T={Tmax},N={NL},hz={hz},p={p},mtg-n=2.dat', E3),
        (f'HtestIX,T={Tmax},N={NL},hz={hz},p={p},mtg-n=3.dat', E4)
    ]

    # 批量写入文件
    for filename, data in files_data:
        with open(os.path.join(base_dir, filename), 'w') as file:
            for i in range(len(data)):
                file.write(f"{deltaT * (i + 1)},{np.real(holevo_information(data[i]))}\n")


def main():
    NL = 7
    J = 1
    hx = 1
    hz = 0.0
    Tmax = 10
    deltaT = 0.1
    p = 0.001

    # 文件路径
    base_dir = "results"
    os.makedirs(base_dir, exist_ok=True)
    num_cores = mp.cpu_count() - 1
    time_start = time.time()
    pool = mp.Pool(processes=num_cores)
    compute_func = partial(compute_time_step, NL=NL, J=J, hx=hx, hz=hz, deltaT=deltaT, p=p)
    tasks = range(int(Tmax / deltaT))
    results = pool.map(compute_func, tasks)
    pool.close()
    pool.join()
    write_results(base_dir, deltaT, results, Tmax, NL, hz, p)
    time_end = time.time()
    print(f"总计算时间: {time_end - time_start:.2f} 秒")

if __name__ == "__main__":
    main()
