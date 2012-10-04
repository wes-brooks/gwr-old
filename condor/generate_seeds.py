import random
N=10000

f = open("seeds.txt", "w")

for i in range(N):
    r = random.random()
    f.write(str(r) + '\n')


