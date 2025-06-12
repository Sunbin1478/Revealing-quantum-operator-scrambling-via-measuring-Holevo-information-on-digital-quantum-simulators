import multiprocessing as mp
from functools import partial
import os
os.environ["QIBO_BACKEND"] = "numpy"  # 强制使用 NumPy 后端，避免使用 CuPy
import numpy as np
import time
from qibo.models import Circuit
from qibo import gates, models
def list(p):
    zip =(["X",p/3],["Y",p/3],["Z",p/3])
    return zip
def bellst(c, n, m, p):
    c.add(gates.H(n))
    c.add(gates.PauliNoiseChannel(n, list(p)))

    c.add(gates.CNOT(n, m))
    c.add(gates.PauliNoiseChannel(n, list(p)))
    c.add(gates.PauliNoiseChannel(m, list(p)))
    return c

def ibellst(c, n, m, p):
    c.add(gates.CNOT(n, m))
    c.add(gates.PauliNoiseChannel(n,list(p)))
    c.add(gates.PauliNoiseChannel(m, list(p)))

    c.add(gates.H(n))
    c.add(gates.PauliNoiseChannel(n, list(p)))

    return c


# 对演化的末态做个操作，使得bell态的信息容易读取，|00>,|10>,|01>,|11>
# ==========================================================================
def RZZ(J, dt, c, q1, q2, p):
    c.add(gates.CNOT(q1, q2))
    c.add(gates.PauliNoiseChannel(q1, list(p)))
    c.add(gates.PauliNoiseChannel(q2, list(p)))

    c.add(gates.RZ(q2, theta=-2 * J * dt))
    c.add(gates.PauliNoiseChannel(q2, list(p)))

    c.add(gates.CNOT(q1, q2))
    c.add(gates.PauliNoiseChannel(q1,list(p)))
    c.add(gates.PauliNoiseChannel(q2, list(p)))
    return c
def Hmix_dt_evl_up(c, NL, J, hz, hx, dt, p):
    for i in range(int(NL / 2)):
        c = RZZ(J, dt, c, 2 * i, 2 * i + 1, p)
    for i in range(int((NL - 1) / 2)):
        c = RZZ(J, dt, c, 2 * i + 1, 2 * i + 2, p)
    for i in range(NL):
        c.add(gates.RZ(i, theta=-2 * hz * dt))
        c.add(gates.PauliNoiseChannel(i, list(p)))
    for i in range(NL):
        c.add(gates.RX(i, theta=-2 * hx * dt))
        c.add(gates.PauliNoiseChannel(i, list(p)))
    return c

def Hmix_dt_evl_down(c, NL, J, hz, hx, dt, p):
    for i in range(int(NL / 2)):
        c = RZZ(J, dt, c, 2 * i, 2 * i + 1, p)
    for i in range(int((NL - 1) / 2)):
        c = RZZ(J, dt, c,  2 * i + 1,  2 * i + 2, p)
    for i in range(NL):
        c.add(gates.RZ( i, theta=-2 * hz * dt))
        c.add(gates.PauliNoiseChannel( i, list(p)))
    for i in range(NL):
        c.add(gates.RX( i, theta=-2 * hx * dt))
        c.add(gates.PauliNoiseChannel( i, list(p)))
    return c
def time_evl_get_Lbar(NL, J, hx, hz, T, ci, p, I, gate_type):
    step = max(1, int(T / 0.1))
    dt = T / step
    p = ci * p
    c = models.Circuit(2*NL, density_matrix=True)
    for i in range(NL):
        c = bellst(c, i, i + NL, p)
    for _ in range(step):
        c = Hmix_dt_evl_down(c, NL, J, hz, hx, dt, p)

    if gate_type == 'X':
        c.add(gates.X(int(0)))
    else:
        c.add(gates.Y(int(0)))


    c.add(gates.PauliNoiseChannel(int(NL / 2), list(p)))

    for _ in range(step):
        c = Hmix_dt_evl_up(c, NL, J, hz, hx, -dt, p)

    c = ibellst(c, I, NL+I, p)
    result = c()
    return oplbar(result, I, NL+I)
def entropy(p):
    return -np.sum(p * np.log2(p + 1e-14))  # 添加小量以避免log(0)
def holevo_information(data1I, data2I, data1X, data2X, data1Y, data2Y, data1Z, data2Z):

    avg_I = (data1I + data2I) / 2
    avg_X = (data1X + data2X) / 2
    avg_Y = (data1Y + data2Y) / 2
    avg_Z = (data1Z + data2Z) / 2

    avg_entropy = entropy(np.array([avg_I, avg_X, avg_Y, avg_Z]))
    entropy1 = entropy(np.array([data1I, data1X, data1Y, data1Z]))
    entropy2 = entropy(np.array([data2I, data2X, data2Y, data2Z]))
    holevo_info = avg_entropy - 0.5 * (entropy1 + entropy2)
    return holevo_info
def oplbar(result, NL, I):
    A = result.state()
    qubit_indices = [I, NL]

    prob_00 = probability_00(A, qubit_indices)
    prob_01 = probability_01(A, qubit_indices)
    prob_10 = probability_10(A, qubit_indices)
    prob_11 = probability_11(A, qubit_indices)

    return prob_00, prob_01, prob_10, prob_11


def probability_00(state, qubit_indices):
    return probability(state, qubit_indices, '00')


def probability_01(state, qubit_indices):
    return probability(state, qubit_indices, '01')


def probability_10(state, qubit_indices):
    return probability(state, qubit_indices, '10')


def probability_11(state, qubit_indices):
    return probability(state, qubit_indices, '11')


def probability(state, qubit_indices, target_state):
    num_qubits = int(np.log2(len(state)))
    total_prob = 0.0

    for i in range(len(state)):
        binary_state = f"{i:0{num_qubits}b}"
        if all(binary_state[qubit_indices[j]] == target_state[j] for j in range(len(qubit_indices))):
            total_prob += np.real(state[i, i])

    return total_prob
def compute_time_step(i, NL, J, hx, hz, deltaT, p):
    T = deltaT * (i + 1)
    I = 0
    start_time = time.time()
    E0X_00, E0X_01, E0X_10, E0X_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 1, 1e-10, I, 'X')
    E1X_00, E1X_01, E1X_10, E1X_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 1, p, I, 'X')
    E2X_00, E2X_01, E2X_10, E2X_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 2, p, I, 'X')
    E3X_00, E3X_01, E3X_10, E3X_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 3, p, I, 'X')
    E4X_00, E4X_01, E4X_10, E4X_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 4, p, I, 'X')
    E0Y_00, E0Y_01, E0Y_10, E0Y_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 1, 1e-10, I, 'Y')
    E1Y_00, E1Y_01, E1Y_10, E1Y_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 1, p, I, 'Y')
    E2Y_00, E2Y_01, E2Y_10, E2Y_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 2, p, I, 'Y')
    E3Y_00, E3Y_01, E3Y_10, E3Y_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 3, p, I, 'Y')
    E4Y_00, E4Y_01, E4Y_10, E4Y_11 = time_evl_get_Lbar(NL, J, hx, hz, T, 4, p, I, 'Y')
    E0_result = [E0X_00, E0Y_00, E0X_01, E0Y_01, E0X_10, E0Y_10, E0X_11, E0Y_11]
    E1_result = [E1X_00, E1Y_00, E1X_01, E1Y_01, E1X_10, E1Y_10, E1X_11, E1Y_11]
    E2_result = [
        2 * E1X_00 - E2X_00, 2 * E1Y_00 - E2Y_00,
        2 * E1X_01 - E2X_01, 2 * E1Y_01 - E2Y_01,
        2 * E1X_10 - E2X_10, 2 * E1Y_10 - E2Y_10,
        2 * E1X_11 - E2X_11, 2 * E1Y_11 - E2Y_11
    ]
    E3_result = [
        3 * E1X_00 - 3 * E2X_00 + E3X_00,
        3 * E1Y_00 - 3 * E2Y_00 + E3Y_00,
        3 * E1X_01 - 3 * E2X_01 + E3X_01,
        3 * E1Y_01 - 3 * E2Y_01 + E3Y_01,
        3 * E1X_10 - 3 * E2X_10 + E3X_10,
        3 * E1Y_10 - 3 * E2Y_10 + E3Y_10,
        3 * E1X_11 - 3 * E2X_11 + E3X_11,
        3 * E1Y_11 - 3 * E2Y_11 + E3Y_11
    ]
    E4_result = [
        4 * E1X_00 - 6 * E2X_00 + 4 * E3X_00 - E4X_00,
        4 * E1Y_00 - 6 * E2Y_00 + 4 * E3Y_00 - E4Y_00,
        4 * E1X_01 - 6 * E2X_01 + 4 * E3X_01 - E4X_01,
        4 * E1Y_01 - 6 * E2Y_01 + 4 * E3Y_01 - E4Y_01,
        4 * E1X_10 - 6 * E2X_10 + 4 * E3X_10 - E4X_10,
        4 * E1Y_10 - 6 * E2Y_10 + 4 * E3Y_10 - E4Y_10,
        4 * E1X_11 - 6 * E2X_11 + 4 * E3X_11 - E4X_11,
        4 * E1Y_11 - 6 * E2Y_11 + 4 * E3Y_11 - E4Y_11
    ]

    end_time = time.time()

    return i, E0_result, E1_result, E2_result, E3_result, E4_result
def write_results(base_dir, deltaT, results, Tmax, NL, hz, p):
    sorted_results = sorted(results)
    E0 = np.array([r[1] for r in sorted_results])
    E1 = np.array([r[2] for r in sorted_results])
    E2 = np.array([r[3] for r in sorted_results])
    E3 = np.array([r[4] for r in sorted_results])
    E4 = np.array([r[5] for r in sorted_results])
    files_data = [
        (f'HXY,T={Tmax},N={NL},hz={hz},p={1e-10},mtg-n=0.dat', E0),
        (f'HXY,T={Tmax},N={NL},hz={hz},p={p},mtg-n=0.dat', E1),
        (f'HXY,T={Tmax},N={NL},hz={hz},p={p},mtg-n=1.dat', E2),
        (f'HXY,T={Tmax},N={NL},hz={hz},p={p},mtg-n=2.dat', E3),
        (f'HXY,T={Tmax},N={NL},hz={hz},p={p},mtg-n=3.dat', E4)
    ]

    for filename, data in files_data:
        with open(os.path.join(base_dir, filename), 'w') as file:
            for i in range(len(data)):
                holevo_value = np.real(holevo_information(*data[i]))
                file.write(f"{deltaT * (i + 1)},{holevo_value}\n")
        print(f"已写入文件: {filename}")
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