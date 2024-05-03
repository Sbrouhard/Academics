
total = 14

def fibbonacci(n):
    if n == 0 or n == 1:
        return 1
    else:
        return fibbonacci(n-2) + fibbonacci(n-1)
    
def fibonnaci_dictionary(distance_out):
    fib_dictionary = {}
    for i in range(0, distance_out):
        fib_dictionary[i] = fibbonacci(i)
    
    return fib_dictionary


def get_divisors(n):
    divisors_list = []
    for i in range(1, n):
        if n % i == 0:
            divisors_list.append(i)
    return divisors_list


def get_total_fixed(p):
    fib_dic = fibonnaci_dictionary(2 * p + 1)
    return abs((fib_dic[2 * p] - 1) * (fib_dic[2 * p - 2] - 1) - (fib_dic[2 * p - 1] * (fib_dic[2 * p - 1])))

num_fixed_points_saved = {}
for i in range(1,total):
    num_fixed_points_saved[i] = get_total_fixed(i)


num_unique_fixed_points_saved = {}
def get_unique_fixed_points(p):
    if p == 1:
        num_unique_fixed_points_saved[1] = 1
        return 1
    try:
        return num_unique_fixed_points_saved[p]
    except:
        total = num_fixed_points_saved[p]
        divisors = get_divisors(p)
        for divisor in divisors:
            total -= num_unique_fixed_points_saved[divisor]
        num_unique_fixed_points_saved[p] = total

        return total



print(num_fixed_points_saved[4])

for i in range(1, total):
    get_unique_fixed_points(i)

for i in range(1, total):
    print(f'{i} &  {num_fixed_points_saved[i]} & {get_divisors(i)} & {num_fixed_points_saved[i] - num_unique_fixed_points_saved[i]} & {num_unique_fixed_points_saved[i]} & {num_unique_fixed_points_saved[i] / i} \\\\')