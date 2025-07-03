import time
import math
import random
import csv
from sympy import isprime


def miller_rabin(n, k=5):
    """בודק האם מספר ראשוני באמצעות מבחן מילר-רבין"""
    if n < 2:
        return False, 0
    if n in (2, 3):
        return True, 1
    if n % 2 == 0:
        return False, 1

    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2

    iterations = 0

    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n)
        iterations += 1

        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            iterations += 1
            if x == n - 1:
                break
        else:
            return False, iterations
    return True, iterations


def pollards_rho_improved(n, max_iter=100000, alternative_c=[1, 2, 3, 5, 7, 11]):
    """פירוק לגורמים באמצעות שיטת רו של פולרד"""
    if isprime(n):
        return "Prime", 0
    if n % 2 == 0:
        return 2, 1

    for c in alternative_c:
        x, y = random.randint(1, n - 1), random.randint(1, n - 1)
        d = 1
        iterations = 0

        while d == 1 and iterations < max_iter:
            x = (x ** 2 + c) % n
            y = (y ** 2 + c) % n
            y = (y ** 2 + c) % n
            d = math.gcd(abs(x - y), n)
            iterations += 1
            if iterations >= max_iter:
                break

        if d != 1 and d != n:
            return d, iterations

    return "Timeout", max_iter


def fermat_factorization_improved(n):
    """פירוק לגורמים באמצעות שיטת פרמה"""
    if isprime(n):
        return "Prime", n, 0
    if n % 2 == 0:
        return 2, n // 2, 1

    x = math.isqrt(n) + 1
    iterations = 0

    while True:
        y2 = x * x - n
        y = int(math.sqrt(y2))
        iterations += 1

        if y * y == y2:
            return x - y, x + y, iterations

        x += 1


def quadratic_sieve(n, max_attempts=100000):
    """פירוק לגורמים באמצעות שיטת הסינון הריבועי"""
    if isprime(n):
        return "Prime", n, 0

    x = math.isqrt(n) + 1
    attempts = 0

    while attempts < max_attempts:
        x2 = x * x - n
        y = int(math.sqrt(x2))
        attempts += 1

        if y * y == x2:
            return x - y, x + y, attempts

        x += 1

    return "Timeout", n, max_attempts


# רשימת מספרים לבדיקה
numbers = [8633, 809009, 92296873, 88169891, 4601, 91, 8051, 7031, 2701,
           9509, 13561, 8777, 14429, 12403, 14527, 10123, 12449, 9353,
           25511, 17873, 15, 21, 401, 4507, 59651, 296627]

# רשימה לאחסון התוצאות
results = []

for num in numbers:
    start_time = time.time()
    is_prime, iterations_mr = miller_rabin(num)
    elapsed_time_mr = time.time() - start_time

    start_time = time.time()
    factor_rho, iterations_rho = pollards_rho_improved(num)
    elapsed_time_rho = time.time() - start_time

    start_time = time.time()
    factors_fermat = fermat_factorization_improved(num)
    elapsed_time_fermat = time.time() - start_time

    start_time = time.time()
    factors_sieve = quadratic_sieve(num)
    elapsed_time_sieve = time.time() - start_time

    print(f"Number: {num}, "
          f"Miller-Rabin: Prime={is_prime}, Iterations={iterations_mr}, Time={elapsed_time_mr:.6f}s, "
          f"Pollard's Rho: Factor={factor_rho}, Iterations={iterations_rho}, Time={elapsed_time_rho:.6f}s, "
          f"Fermat: Factors={factors_fermat[0]}, {factors_fermat[1]}, Iterations={factors_fermat[2]}, Time={elapsed_time_fermat:.6f}s, "
          f"Quadratic Sieve: Factors={factors_sieve[0]}, {factors_sieve[1]}, Iterations={factors_sieve[2]}, Time={elapsed_time_sieve:.6f}s")

    results.append([
        num,
        "Prime" if is_prime else "Composite", iterations_mr, elapsed_time_mr,
        factor_rho, iterations_rho, elapsed_time_rho,
        factors_fermat[0], factors_fermat[1], factors_fermat[2], elapsed_time_fermat,
        factors_sieve[0], factors_sieve[1], factors_sieve[2], elapsed_time_sieve
    ])

# שמירת התוצאות לקובץ CSV
csv_filename = "factorization_results.csv"
with open(csv_filename, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow([
        "Number", "Miller-Rabin Prime", "Miller-Rabin Iterations", "Miller-Rabin Time",
        "Pollard's Rho Factor", "Pollard's Rho Iterations", "Pollard's Rho Time",
        "Fermat Factor 1", "Fermat Factor 2", "Fermat Iterations", "Fermat Time",
        "Quadratic Sieve Factor 1", "Quadratic Sieve Factor 2", "Quadratic Sieve Iterations", "Quadratic Sieve Time"
    ])
    writer.writerows(results)

print(f"Results saved to {csv_filename}")
