
GAMMA = 2.08597893692734
K_MIN = 2
K_MAX = 100
sum_x = 0
for i in range(K_MIN, K_MAX):
    sum_x = sum_x + i**(-GAMMA)
NORMALIZING_FACTOR = sum_x
print("NORMALIZING_FACTOR=", NORMALIZING_FACTOR)