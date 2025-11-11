from itertools import product
from sympy import factorint, isprime, gcd
from sage.all import euler_phi, factor
import time

def is_prime(n: int) -> bool:
    if n < 2 or (n % 2 == 0 and n != 2):
        return False
    
    for i in range(3, int(n ** 0.5), 2):
        if n % i == 0:
            return False
    return True

def eratosthenes_sieve(n: int) -> list[int]:
    arr = [i for i in range(n + 1)]
    arr[1] = 0

    for i in range(2, int(n ** 0.5)):
        if is_prime(i):
            for j in range(i ** 2, n + 1, i):
                arr[j] = 0
    arr_of_prime = [i for i in arr if i != 0]

    return arr_of_prime

list_of_primes = eratosthenes_sieve(10 ** 6)

def palindromic_squares_and_circular_primes() -> tuple[list[int], list[int]]:
    arr_of_palindroms = []
    arr_of_primes = []

    for a in range(10 ** 5):
        if a == int(str(a)[::-1]):
            if a ** 2 == int(str(a ** 2)[::-1]):
                arr_of_palindroms.append(a)
    
    digits_can_broke = ["0", "2", "4", "6", "8", "5"] #Числа, наличие которых в нашем числе точно не позволит выполниться условию
    for p in list_of_primes:
        for j in str(p):
            if j in digits_can_broke:
                break
        else:
            cycle_p = str(p)
            while cycle_p != p:
                cycle_p = int(str(cycle_p)[1:] + str(cycle_p)[0])
                if not is_prime(cycle_p):
                    break
                else:
                    arr_of_primes.append(p)
    return (arr_of_palindroms, arr_of_primes)

def palindromic_cubes_and_palindromic_primes() -> tuple[list[int], list[int]]:
    arr_of_palindroms = []
    arr_of_prime_palindroms = []

    for a in range(10 ** 5):
        if a == int(str(a)[::-1]):
            if a ** 3 == int(str(a ** 3)[::-1]):
                arr_of_palindroms.append(a)
    i = 0
    while list_of_primes[i] < 10 ** 5:
        p = list_of_primes[i]
        i += 1
        if p == int(str(p)[::-1]):
            arr_of_prime_palindroms.append(p)
    return (arr_of_palindroms, arr_of_prime_palindroms)

def primes_with_two_digits() -> dict[str, list[int]]:
    dict_of_primes = {}
    dict_of_primes["13"] = []
    dict_of_primes["15"] = []
    dict_of_primes["17"] = []
    dict_of_primes["19"] = []

    lenght = 1
    while len(dict_of_primes["13"]) + len(dict_of_primes["15"]) + len(dict_of_primes["17"]) + len(dict_of_primes["19"]) < 400:
        for digits in [["1", "3"], ["1", "5"], ["1", "7"], ["1", "9"]]:
            if len(dict_of_primes["1" + digits[1]]) < 100:
                for p in product(digits, repeat=lenght):
                    number = int("".join(p))
                    if number % 10 != 5:
                        if isprime(int(number)):
                            dict_of_primes["1" + digits[1]].append(number)
        lenght += 1
    return dict_of_primes

def factorial_plus_one_factors() -> dict[int, dict[int, int]]:
    result_dict = dict()
    number = 1
    max_dig = 0
    arr_of_big_dig = []

    for i in range(2, 51):
        number *= i
        number += 1
        facs = factor(number)
        result_dict[i] = facs
        number -= 1

        max_dig = max(len(result_dict[i]), max_dig)

        if max(result_dict[i][0]) > 10 ** 6:
            arr_of_big_dig.append(number + 1)

    print(max_dig, arr_of_big_dig)

    return result_dict

def euler_phi_direct(n: int) -> int:
    """Вычисляет (n) прямым перебором."""

    Zm_star = [i for i in range(1, n) if gcd(i, n)]

    return len(Zm_star)

def euler_phi_factor(n: int) -> int:
    """Вычисляет (n) через разложение на простые множители."""

    factorization = factorint(n)
    phi = n

    for i in factorization:
        phi *= (1 - 1 / i)

    return int(phi)

def compare_euler_phi_methods(test_values: list[int]) -> dict:
    """
    Сравнивает время работы трёх методов на заданных значениях.
    Возвращает словарь с тремя списками времён (в секундах).
    """
    time_manager = dict()

    start_direct = time.perf_counter()
    euler_phi_direct(test_values)
    end_direct = time.perf_counter()
    time_manager["Замер прямого подсчета"] = end_direct - start_direct

    start_factor = time.perf_counter()
    euler_phi_factor(test_values)
    end_factor = time.perf_counter()
    time_manager["Замер подсчета через факторизацию"] = end_factor - start_factor

    start_sage = time.perf_counter()
    euler_phi(test_values)
    end_sage = time.perf_counter()
    time_manager["Замер счета через sage"] = end_sage - start_sage

    return time_manager

#palindromic_cubes_and_palindromic_primes()
#palindromic_squares_and_circular_primes()
#print(primes_with_two_digits())
#print(compare_euler_phi_methods(10000324))
print(factorial_plus_one_factors())
print(euler_phi_factor(134517 * 7))
print(euler_phi(134517 * 7))
