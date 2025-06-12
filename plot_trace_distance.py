import random
import numpy as np
from numpy import sqrt
from numpy.linalg import qr
from projectq import MainEngine
from projectq.ops import BasicGate, All, Measure, Rx, Ry, T, CZ
from projectq.backends import Simulator
import matplotlib.pyplot as plt
from scipy.linalg import sqrtm
import json

# 设置全局绘图字体为 Times New Roman
plt.rcParams.update({
    'font.size': 40,  # 全局字体大小
    'axes.titlesize': 36,  # 标题字体大小
    'axes.labelsize': 36,  # 坐标轴标签字体大小
    'legend.fontsize': 28,  # 图例字体大小
    'xtick.labelsize': 36,  # X 轴刻度标签字体大小
    'ytick.labelsize': 36,  # Y 轴刻度标签字体大小
    'font.family': 'Times New Roman'  # 设置字体族为 Times New Roman
})


# 从 JSON 文件加载门记录的函数
def load_gate_records(L):
    filename = f'gate_records_N={L + 1}_hz=1e-05.json'
    try:
        with open(filename, 'r') as file:
            gate_records = json.load(file)
        return gate_records
    except FileNotFoundError:
        print(f"文件 {filename} 未找到。")
        return []


# 初始化量子计算引擎
eng = MainEngine(backend=Simulator())
colors = plt.cm.tab20(np.linspace(0, 1, 15))  # 使用原始颜色方案
markers = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x', 'd', '|']
plt.figure(figsize=(12, 9))


# 生成最大混合态密度矩阵的函数
def Maximum_mixed_state(L):
    Maximum_mixed_density = np.eye(2 ** L) / (2 ** L)
    return Maximum_mixed_density


# 应用单个门记录到量子位的函数
def apply_gate_record(qubits, gate_record):
    cycle_number, qubit_index, gate_name, angle = gate_record
    if gate_name == 'Rx':
        Rx(angle) | qubits[qubit_index]
    elif gate_name == 'Ry':
        Ry(angle) | qubits[qubit_index]
    elif gate_name == 'T':
        T | qubits[qubit_index]
    elif gate_name == 'CZ':
        CZ | (qubits[qubit_index], qubits[qubit_index + 1])

def apply_gate_records(qubits, gate_records):
    for record in gate_records:
        apply_gate_record(qubits, record)


num_cycles = 100

def random_unitarygate(L, gate_records):
    qubits = eng.allocate_qureg(L)
    apply_gate_records(qubits, gate_records)
    return qubits

def hilbert_schmidt_distance(rho, sigma):
    difference = rho - sigma
    return np.linalg.norm(difference, 'fro')


for L in range(3, 7, 1):
    gate_records_all = load_gate_records(L)
    if not gate_records_all:
        continue
    trace_distances = []
    sample_counts = []
    data = np.zeros((2 ** L, 2 ** L), dtype=np.complex128)

    for m in range(1, 1001, 1):
        # 为每个样本初始化量子态
        eng.flush()
        qubits = random_unitarygate(L, gate_records_all[m % len(gate_records_all)])
        eng.flush()
        indices, cheat_state = eng.backend.cheat()  # 获取状态向量
        cheat_state_array = np.array(cheat_state)
        cheat_state_column = cheat_state_array[:, np.newaxis]
        cheat_state_dagger = np.conj(cheat_state_column).T
        cheat_density = np.dot(cheat_state_column, cheat_state_dagger)  # 纯态密度矩阵
        data += cheat_density  # 累加密度矩阵
        All(Measure) | qubits  # 清除量子位以进行下一次迭代
        data1 = data / m
        Maximum_mixed_density = Maximum_mixed_state(L)  # 最大混合态矩阵
        a = hilbert_schmidt_distance(data1, Maximum_mixed_density)  # 计算迹距离
        trace_distances.append(a.real)
        sample_counts.append(m)
    ln_trace_distances = np.log10(trace_distances)
    ln_sample_counts = np.log10(sample_counts)
    plt.scatter(
        ln_sample_counts,
        ln_trace_distances,
        label=f'L={L}',
        color=colors[L - 1],
        marker=markers[L - 1],
        s=50
    )
# 绘制理论曲线
sample_counts = np.array(sample_counts)
ln_sample_counts1 = np.log10(sample_counts)
ln_theoretical_trace_distances = -0.5 * ln_sample_counts1
plt.plot(
    ln_sample_counts,
    ln_theoretical_trace_distances,
    label='Theoretical Curve',
    color='black',
    linestyle='--',
    linewidth=2
)

plt.xlabel('log$_{10}$(M)', fontsize=32)
plt.ylabel(r'$\log_{10}D(\rho^{\mathrm{est}},\ \rho^{\mathrm{max}})$', fontsize=32)
plt.title(
    r'$\log_{10}D(\rho^{\mathrm{est}},\ \rho^{\mathrm{max}})$ vs. log$_{10}$(M) for Different L',
    fontsize=32
)
plt.xticks(
    ticks=np.linspace(ln_sample_counts.min(), ln_sample_counts.max(), 4)
)

plt.yticks(ticks=[-0.0, -0.5, -1.0, -1.5])

plt.legend(fontsize=28, loc='upper right')
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('enhanced_plot_times_new_roman.png', dpi=300, bbox_inches='tight')
plt.show()



