import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 定义标准化边的函数，确保边的表示唯一
def normalize_edge(p1, p2):
    return tuple(sorted([p1, p2]))

# 获取路径上相邻立方体共面上的四条棱
def get_shared_face_edges(c1, c2):
    i1, j1, k1 = c1
    i2, j2, k2 = c2
    di, dj, dk = i2 - i1, j2 - j1, k2 - k1
    
    if abs(di) + abs(dj) + abs(dk) != 1:
        raise ValueError("立方体不相邻")
    
    if di == 1:
        x = i1 + 1
        edges = [
            normalize_edge((x, j1, k1), (x, j1+1, k1)),
            normalize_edge((x, j1, k1), (x, j1, k1+1)),
            normalize_edge((x, j1+1, k1), (x, j1+1, k1+1)),
            normalize_edge((x, j1, k1+1), (x, j1+1, k1+1)),
        ]
    elif di == -1:
        x = i1
        edges = [
            normalize_edge((x, j1, k1), (x, j1+1, k1)),
            normalize_edge((x, j1, k1), (x, j1, k1+1)),
            normalize_edge((x, j1+1, k1), (x, j1+1, k1+1)),
            normalize_edge((x, j1, k1+1), (x, j1+1, k1+1)),
        ]
    elif dj == 1:
        y = j1 + 1
        edges = [
            normalize_edge((i1, y, k1), (i1+1, y, k1)),
            normalize_edge((i1, y, k1), (i1, y, k1+1)),
            normalize_edge((i1+1, y, k1), (i1+1, y, k1+1)),
            normalize_edge((i1, y, k1+1), (i1+1, y, k1+1)),
        ]
    elif dj == -1:
        y = j1
        edges = [
            normalize_edge((i1, y, k1), (i1+1, y, k1)),
            normalize_edge((i1, y, k1), (i1, y, k1+1)),
            normalize_edge((i1+1, y, k1), (i1+1, y, k1+1)),
            normalize_edge((i1, y, k1+1), (i1+1, y, k1+1)),
        ]
    elif dk == 1:
        z = k1 + 1
        edges = [
            normalize_edge((i1, j1, z), (i1+1, j1, z)),
            normalize_edge((i1, j1, z), (i1, j1+1, z)),
            normalize_edge((i1+1, j1, z), (i1+1, j1+1, z)),
            normalize_edge((i1, j1+1, z), (i1+1, j1+1, z)),
        ]
    elif dk == -1:
        z = k1
        edges = [
            normalize_edge((i1, j1, z), (i1+1, j1, z)),
            normalize_edge((i1, j1, z), (i1, j1+1, z)),
            normalize_edge((i1+1, j1, z), (i1+1, j1+1, z)),
            normalize_edge((i1, j1+1, z), (i1+1, j1+1, z)),
        ]
    return edges

# 示例路径（假设有三条路径）
paths = [
    [(0,0,0), (1,0,0), (2,0,0)],  # 第一条路径
    [(0,1,0), (0,2,0), (0,3,0)],  # 第二条路径
    [(0,0,1), (1,0,1), (2,0,1)]   # 第三条路径
]
path_colors = ['r', 'g', 'b']  # 路径颜色

# 创建3D图形
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 生成所有可能的边（6x6x6立方体网格，顶点从0到6）
all_edges = set()
# x方向的边
for i in range(6):
    for j in range(7):
        for k in range(7):
            p1 = (i, j, k)
            p2 = (i+1, j, k)
            all_edges.add(normalize_edge(p1, p2))
# y方向的边
for i in range(7):
    for j in range(6):
        for k in range(7):
            p1 = (i, j, k)
            p2 = (i, j+1, k)
            all_edges.add(normalize_edge(p1, p2))
# z方向的边
for i in range(7):
    for j in range(7):
        for k in range(6):
            p1 = (i, j, k)
            p2 = (i, j, k+1)
            all_edges.add(normalize_edge(p1, p2))

# 收集路径上需要省略的边
omit_edges = set()
for path in paths:
    for idx in range(len(path) - 1):
        c1 = path[idx]
        c2 = path[idx + 1]
        edges = get_shared_face_edges(c1, c2)
        omit_edges.update(edges)

# 需要绘制的边
draw_edges = all_edges - omit_edges

# 绘制迷宫的边缘
for edge in draw_edges:
    p1, p2 = edge
    ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='black', linewidth=1)

# 绘制路径
for path, color in zip(paths, path_colors):
    path_centers = [(i+0.5, j+0.5, k+0.5) for i, j, k in path]
    ax.plot([p[0] for p in path_centers], [p[1] for p in path_centers], [p[2] for p in path_centers], 
            color=color, linewidth=2)

# 设置图形参数
ax.view_init(elev=30, azim=45)
ax.set_xlim(0, 6)
ax.set_ylim(0, 6)
ax.set_zlim(0, 6)
ax.set_axis_off()
plt.show()