import matplotlib.pyplot as plt
import numpy as np
# 设置 Times New Roman 字体
plt.rcParams['font.family'] = 'Times New Roman'
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

original_cmap = plt.cm.viridis

colors = original_cmap(np.linspace(0, 1, 256))
colors[-1, :3] = colors[-1, :3] * 1
custom_cmap = LinearSegmentedColormap.from_list('viridis_modified', colors)
# 读取文件内容的函数
def filedata(filename):
    times = []
    grid_point_1 = []
    grid_point_2 = []
    grid_point_3 = []
    grid_point_4 = []
    grid_point_5 = []
    grid_point_6 = []
    grid_point_7 = []
    with open(filename, "r") as file:
        lines = file.readlines()
        for i in range(0, len(lines), 9):  # 每个时间步有6行数据和一个空行
            if i + 7 >= len(lines):
                break  # 防止索引越界
            time_line = lines[i].strip()
            grid_point_1_line = lines[i + 1].strip()
            grid_point_2_line = lines[i + 2].strip()
            grid_point_3_line = lines[i + 3].strip()
            grid_point_4_line = lines[i + 4].strip()
            grid_point_5_line = lines[i + 5].strip()
            grid_point_6_line = lines[i + 6].strip()
            grid_point_7_line = lines[i + 7].strip()
            # 提取时间和格点值
            time_str = time_line.split()[1][:-1]
            times.append(float(time_str))
            grid_point_1_value = float(grid_point_1_line.split(":")[1])
            grid_point_2_value = float(grid_point_2_line.split(":")[1])
            grid_point_3_value = float(grid_point_3_line.split(":")[1])
            grid_point_4_value = float(grid_point_4_line.split(":")[1])
            grid_point_5_value = float(grid_point_5_line.split(":")[1])
            grid_point_6_value = float(grid_point_6_line.split(":")[1])
            grid_point_7_value = float(grid_point_7_line.split(":")[1])

            grid_point_1.append(grid_point_1_value)
            grid_point_2.append(grid_point_2_value)
            grid_point_3.append(grid_point_3_value)
            grid_point_4.append(grid_point_4_value)
            grid_point_5.append(grid_point_5_value)
            grid_point_6.append(grid_point_6_value)
            grid_point_7.append(grid_point_7_value)
    return grid_point_1, grid_point_2, grid_point_3, grid_point_4, grid_point_5,grid_point_6,grid_point_7,times
filename = "holevoIX_information_N=7_hz=0.3.txt"
filename2 = "holevoIX_information_N=7_hz=0.0txt"

grid_point_1, grid_point_2, grid_point_3, grid_point_4, grid_point_5,grid_point_6,grid_point_7,times = filedata(filename)
grid_point_11, grid_point_21, grid_point_31, grid_point_41, grid_point_51,grid_point_61,grid_point_71,times = filedata(filename2)
data_matrix = np.array(
    [grid_point_1, grid_point_2, grid_point_3, grid_point_4, grid_point_5, grid_point_6, grid_point_7]).T
data_matrix1 = np.array(
    [grid_point_11, grid_point_21, grid_point_31, grid_point_41, grid_point_51, grid_point_61, grid_point_71]).T
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 12))  # 增大图形大小以适应更大的字体
vmin = min(data_matrix.min(), data_matrix1.min())
vmax = max(data_matrix.max(), data_matrix1.max())
from matplotlib.colors import PowerNorm
vmin = min(data_matrix.min(), data_matrix1.min())
vmax = max(data_matrix.max(), data_matrix1.max())
norm = PowerNorm(gamma=0.9, vmin=vmin, vmax=1)

# 绘制热图
for i, (ax, data, title) in enumerate(zip([ax1, ax2], [data_matrix1, data_matrix], ['hz=0.0', 'hz=0.3'])):
    cax = ax.imshow(data, aspect='auto', cmap='viridis', origin='lower', norm=norm)

    # 设置横坐标标签
    ax.set_xlabel(r'Site $i$', fontsize=48, labelpad=15) 
    ax.set_xticks([0, 2, 4, 6])
    ax.set_xticklabels(["1", "3", "5", "7"], fontsize=48)
    if i == 0:
        num_y_ticks = 4
        ytick_positions = np.linspace(0, len(times) - 1, num_y_ticks).astype(int)
        ytick_labels = [f"{times[j]:.1f}" for j in ytick_positions]
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels(ytick_labels, fontsize=48)
        ax.set_ylabel('Time', fontsize=54, labelpad=5)
    else:
        ax.set_yticks([])

    ax1.set_title('(a) hz=0.0', fontsize=54, pad=20, loc='center')
    ax2.set_title('(b) hz=0.3', fontsize=54, pad=20, loc='center')

# 添加统一的colorbar
cbar_ax = fig.add_axes([0.85, 0.25, 0.02, 0.60])
cbar = fig.colorbar(cax, cax=cbar_ax, ticks=[0, 1])
cbar.ax.tick_params(labelsize=48)
cbar.set_label('Holevo Information ($\chi$)', fontsize=48, labelpad=5)
plt.tight_layout(rect=[0.05, 0.10, 0.80, 0.95])
plt.show()




