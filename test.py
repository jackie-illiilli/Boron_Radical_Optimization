import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import deque

# 定义迷宫大小
N = 6
np.random.seed(2)
# 定义墙壁颜色和透明度
wall_color = 'gray'
wall_alpha = 0.1

# 定义路径颜色
path_colors = ['red', 'green', 'blue']

# 定义入口和出口
starts = [(0,0,0), (0,2,2), (0,5,5)]
ends = [(5,0,0), (5,2,2), (5,5,5)]
# starts = [(0,2,2)]
# ends = [(5,2,2)]
# path_colors = ['green']
starts_set = set(starts)
ends_set = set(ends)

# 初始化访问标记和墙壁数组
visited = np.zeros((N, N, N), dtype=bool)
walls_x = np.ones((N-1, N, N), dtype=bool)  # x方向墙壁，介于i和i+1之间
walls_y = np.ones((N, N-1, N), dtype=bool)  # y方向墙壁，介于j和j+1之间
walls_z = np.ones((N, N, N-1), dtype=bool)  # z方向墙壁，介于k和k+1之间

# 递归回溯算法生成迷宫
def carve_maze(i, j, k):
    visited[i, j, k] = True
    directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    np.random.shuffle(directions)
    for di, dj, dk in directions:
        ni, nj, nk = i + di, j + dj, k + dk
        if 0 <= ni < N and 0 <= nj < N and 0 <= nk < N and not visited[ni, nj, nk]:
            if di == 1:
                walls_x[i, j, k] = False
            elif di == -1:
                walls_x[ni, j, k] = False
            elif dj == 1:
                walls_y[i, j, k] = False
            elif dj == -1:
                walls_y[i, nj, k] = False
            elif dk == 1:
                walls_z[i, j, k] = False
            elif dk == -1:
                walls_z[i, j, nk] = False
            carve_maze(ni, nj, nk)

# 从第一个入口开始生成迷宫
carve_maze(*starts[0])

# 获取相邻单元格（无墙）
def get_neighbors(i, j, k):
    neighbors = []
    if i > 0 and not walls_x[i-1, j, k]:
        neighbors.append((i-1, j, k))
    if i < N-1 and not walls_x[i, j, k]:
        neighbors.append((i+1, j, k))
    if j > 0 and not walls_y[i, j-1, k]:
        neighbors.append((i, j-1, k))
    if j < N-1 and not walls_y[i, j, k]:
        neighbors.append((i, j+1, k))
    if k > 0 and not walls_z[i, j, k-1]:
        neighbors.append((i, j, k-1))
    if k < N-1 and not walls_z[i, j, k]:
        neighbors.append((i, j, k+1))
    return neighbors

# 使用BFS找到从start到end的路径
def find_path(start, end):
    queue = deque([start])
    parent = {start: None}
    while queue:
        current = queue.popleft()
        if current == end:
            break
        for neighbor in get_neighbors(*current):
            if neighbor not in parent:
                parent[neighbor] = current
                queue.append(neighbor)
    path = []
    current = end
    while current is not None:
        path.append(current)
        current = parent.get(current)
    path.reverse()
    return path

# 找到三条路径
paths = []
for start, end in zip(starts, ends):
    path = find_path(start, end)
    paths.append(path)

# 可视化迷宫
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制内部墙壁
for i in range(N-1):
    for j in range(N):
        for k in range(N):
            if walls_x[i, j, k]:
                yy, zz = np.meshgrid(np.linspace(j, j+1, 2), np.linspace(k, k+1, 2))
                xx = np.full_like(yy, i+1)
                ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

for i in range(N):
    for j in range(N-1):
        for k in range(N):
            if walls_y[i, j, k]:
                xx, zz = np.meshgrid(np.linspace(i, i+1, 2), np.linspace(k, k+1, 2))
                yy = np.full_like(xx, j+1)
                ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

for i in range(N):
    for j in range(N):
        for k in range(N-1):
            if walls_z[i, j, k]:
                xx, yy = np.meshgrid(np.linspace(i, i+1, 2), np.linspace(j, j+1, 2))
                zz = np.full_like(xx, k+1)
                ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

# 绘制外部墙壁
# x=0
for j in range(N):
    for k in range(N):
        if (0,j,k) not in starts_set:
            yy, zz = np.meshgrid(np.linspace(j, j+1, 2), np.linspace(k, k+1, 2))
            xx = np.full_like(yy, 0)
            ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

# x=6
for j in range(N):
    for k in range(N):
        if (5,j,k) not in ends_set:
            yy, zz = np.meshgrid(np.linspace(j, j+1, 2), np.linspace(k, k+1, 2))
            xx = np.full_like(yy, 6)
            ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

# y=0 and y=6
for i in range(N):
    for k in range(N):
        # y=0
        xx, zz = np.meshgrid(np.linspace(i, i+1, 2), np.linspace(k, k+1, 2))
        yy = np.full_like(xx, 0)
        ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)
        # y=6
        yy = np.full_like(xx, 6)
        ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

# z=0 and z=6
for i in range(N):
    for j in range(N):
        # z=0
        xx, yy = np.meshgrid(np.linspace(i, i+1, 2), np.linspace(j, j+1, 2))
        zz = np.full_like(xx, 0)
        ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)
        # z=6
        zz = np.full_like(xx, 6)
        ax.plot_surface(xx, yy, zz, color=wall_color, alpha=wall_alpha)

# 绘制路径
for path, color in zip(paths, path_colors):
    path_centers = [(i+0.5, j+0.5, k+0.5) for i, j, k in path]
    ax.plot([p[0] for p in path_centers], [p[1] for p in path_centers], [p[2] for p in path_centers], color=color, linewidth=5, alpha=0.7)

# 设置视角
ax.view_init(elev=30, azim=45)

# 设置坐标轴范围并隐藏轴
ax.set_xlim(0, N)
ax.set_ylim(0, N)
ax.set_zlim(0, N)
ax.set_axis_off()

plt.show()