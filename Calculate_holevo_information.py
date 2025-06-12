import numpy as np
def entropy(p):
    return -np.sum(p * np.log2(p + 1e-14))  # 添加小量以避免log(0)
def holevo_information(data1I, data2I, data1X, data2X, data1Y, data2Y, data1Z, data2Z):
    arrays = [data1I, data2I, data1X, data2X, data1Y, data2Y, data1Z, data2Z]
    arrays = [np.array(arr) for arr in arrays]
    assert all(arr.shape == (141, 7) for arr in arrays), "All input arrays must be 21x5 in shape"
    data1I, data2I, data1X, data2X, data1Y, data2Y, data1Z, data2Z = arrays
    holevo_info = np.zeros((141, 7))

    for t in range(141):  # 对每个时间点
        for g in range(7):  # 对每个格点
            # 计算平均分布
            avg_I = (data1I[t, g] + data2I[t, g]) / 2
            avg_X = (data1X[t, g] + data2X[t, g]) / 2
            avg_Y = (data1Y[t, g] + data2Y[t, g]) / 2
            avg_Z = (data1Z[t, g] + data2Z[t, g]) / 2
            # 计算平均分布的熵
            avg_entropy = entropy(np.array([avg_I, avg_X, avg_Y, avg_Z]))
            # 计算各个分布的熵
            entropy1 = entropy(np.array([data1I[t, g], data1X[t, g], data1Y[t, g], data1Z[t, g]]))
            entropy2 = entropy(np.array([data2I[t, g], data2X[t, g], data2Y[t, g], data2Z[t, g]]))
            # 计算 Holevo inforamtion
            holevo_info[t, g] = avg_entropy - 0.5 * (entropy1 + entropy2)

    return holevo_info

data1I ="""
"""
data2I = """
"""
data1X = """
"""
data2X = """
"""
data1Y = """
"""
data2Y = """
"""
data1Z = """
"""
data2Z = """
"""


def convert_data(raw_data):
    lines = raw_data.strip().split('\n')
    converted_data = []
    for line in lines:
        line = line.replace('[', '').replace(']', '')
        values = line.split(',')
        row = [float(value.strip()) for value in values if value.strip()]
        if row:
            converted_data.append(row)
    return converted_data

data1I = convert_data(data1I)
data2I = convert_data(data2I)
data1X = convert_data(data1X)
data2X = convert_data(data2X)
data1Y = convert_data(data1Y)
data2Y = convert_data(data2Y)
data1Z = convert_data(data1Z)
data2Z = convert_data(data2Z)
# 计算每个格点在每个时间点上的Holevo信息
holevo_info = holevo_information(data1I, data2I, data1X, data2X, data1Y, data2Y, data1Z, data2Z)
with open("holevoXY_information_N=7_hz=0.0.txt", "w") as f:
    for t, holevo_t in enumerate(holevo_info):
        f.write(f"Time {t * 0.1}:\n")
        for g, holevo_g in enumerate(holevo_t):
            f.write(f"Grid point {g}: {holevo_g}\n")
        f.write("\n")

print("Holevo information has been calculated and stored in the holevoXY_information.txt file")
