import matplotlib.pyplot as plt
for step in range(50):
    x = []
    y = []
    with open('chug' + str(step)) as f:
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
    plt.savefig('foo' + str(step) + '.png')
    plt.clf()
