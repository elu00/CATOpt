
import matplotlib.pyplot as plt
x = []
y = []
with open('chug4999') as f:
    n, m = map(int,f.readline().split())
    for _ in range(n):
        line = f.readline()
        a, b = map(float, line.split())
        x.append(a)
        y.append(b)
    for _ in range(m):
        line = f.readline()
        i, j = map(int, line.split())
        plt.plot([x[i], x[j]], [y[i], y[j]]) 
plt.plot(x, y, 'o')
plt.savefig("final.png")
plt.clf()
