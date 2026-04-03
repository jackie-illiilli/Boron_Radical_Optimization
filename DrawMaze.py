import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import deque

# Define maze size
N = 6
np.random.seed(2)
# Define wall color and transparency
wall_color = 'gray'
wall_alpha = 0.1

# Define path colors
path_colors = ['red', 'green', 'blue']

# Define entrances and exits
starts = [(0,0,0), (0,2,2), (0,5,5)]
ends = [(5,0,0), (5,2,2), (5,5,5)]
# starts = [(0,2,2)]
# ends = [(5,2,2)]
# path_colors = ['green']
starts_set = set(starts)
ends_set = set(ends)

# Initialize visited markers and wall arrays
visited = np.zeros((N, N, N), dtype=bool)
walls_x = np.ones((N-1, N, N), dtype=bool)  # x-direction walls, between i and i+1
walls_y = np.ones((N, N-1, N), dtype=bool)  # y-direction walls, between j and j+1
walls_z = np.ones((N, N, N-1), dtype=bool)  # z-direction walls, between k and k+1

# Recursive backtracking algorithm to generate the maze
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

# Start generating maze from the first entrance
carve_maze(*starts[0])

# Get adjacent cells (no walls)
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

# Use BFS to find the path from start to end
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

# Find three paths
paths = []
for start, end in zip(starts, ends):
    path = find_path(start, end)
    paths.append(path)

# Visualize the maze
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Draw inner walls
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

# Draw outer walls
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

# Draw paths
for path, color in zip(paths, path_colors):
    path_centers = [(i+0.5, j+0.5, k+0.5) for i, j, k in path]
    ax.plot([p[0] for p in path_centers], [p[1] for p in path_centers], [p[2] for p in path_centers], color=color, linewidth=5, alpha=0.7)

# Set viewpoint
ax.view_init(elev=30, azim=45)

# Set axis limits and hide axes
ax.set_xlim(0, N)
ax.set_ylim(0, N)
ax.set_zlim(0, N)
ax.set_axis_off()

plt.show()