from math import gcd
from sage.all import *
from random import randint
from itertools import product, combinations

ISU = 502463
N = ISU % 20

def subgroups_of_Sm(N: int) -> dict:
    m = 4 + N % 5
    
    Sm = SymmetricGroup(m)
    subgroup = Sm.subgroups()

    print(len(subgroup))
    print(list(choice(subgroup)))

    sub_ord = len(Sm) / (N % len(subgroup))
    our_group = 0

    for i in subgroup:
        if len(i) == sub_ord:
            our_group = i

    if our_group == 0:
        return "Подгруппы с нужным индексом нет"
    
    left_classes = [[str(g) for g in coset] for coset in Sm.cosets(our_group, side='left')]
    right_classes = [[str(g) for g in coset] for coset in Sm.cosets(our_group, side='right')]
    
    is_normal = our_group.is_normal(Sm)

    return {
        "Левые классы": left_classes,
        "Правые классы": right_classes,
        "Нормальная ли": is_normal
    }

def element_powers_in_Sm(N: int) -> dict:
    """
    Перебираем элементы группы до тех пор, пока не найдем элемент нужного порядка
    """
    m = 4 + N % 5
    n1, n2, n3 = N % 6, (N + 1) % 6, (N + 2) % 6

    Sm = SymmetricGroup(m)
    result = dict()

    g = 0
    for i in Sm:
        ord = i.order()
        if ord == N:
            g = i
            break
    
    if g == 0:
        print("Вгруппе нет элемента индекса N")
        return 0

    result["gn1"] = (g ** n1).order()
    result["gn2"] = (g ** n2).order()
    result["gn3"] = (g ** n3).order()

    print(g)

    return result

def solve_sigma_power_eq(N: int) -> dict:
    """
    Перебираем элементы группы, пока не найдем решение
    """
    m = 4 + N % 5
    n = 2 + N % 10

    Sm = SymmetricGroup(m)
    
    cycle = "(" + ",".join(str(i) for i in range(1, m)) + ")"

    sigma = Sm(cycle)
    solves = []

    for g in Sm:
        if g ** n == sigma:
            solves.append(g)
    
    random_solutions = []
    for i in range(3):
        random_solutions.append(choice(solves))

    return {
        "solves": solves,
        "rand3": random_solutions
    }

def elements_of_order_k_in_cyclic_group(N: int) -> dict:
    """
    Группа Zm является циклической, порожденной 1, и ее порядок равен m
    """
    n = 2 + N % 10
    m = 4 + N % 10
    k = 1 + N % 7
    
    Zm = Zmod(n)

    ord_k_elem = []
    gk_elem = []

    for g in Zm:
        if g * k % m == 0:
            gk_elem.append(g)
        if g.order() == k:
            ord_k_elem.append(g)

    return {
        "gk=e": gk_elem,
        "ord_k": ord_k_elem
    }

def subgroups_of_Zm_star(N: int) -> list:
    """
    Строим группу Zm и среди всех комбинаций элементов группы находим такие, которые состявляют группу
    """

    m = 4 + N % 5
    #Zm = Integers(m).unit_group()
    Zm = Zmod(m)
    units = [a for a in Zm if gcd(int(a), m) == 1]
    Zm_ord = len(Zm)

    subgroups = [[1]]

    for r in range(2, Zm_ord + 1):
        for i in combinations(units, r):
            sub = list(i)

            if 1 not in sub:
                continue
            
            is_group = True
            for a in sub:
                if a.inverse() not in sub:
                    is_group = False
                    break
                for b in sub:
                    if (a * b) not in sub:
                        is_group = False
                        break
                if not is_group:
                    break
            else:
                subgroups.append(sub)

    return subgroups

def order_of_sr(N: int) -> int:
    p, s, r = 23, 17, 45

    Zp = [i for i in range(1, p) if gcd(i, p) == 1]

    sr = s ** r % p
    cycle_elem = sr
    ord = 1

    while cycle_elem != 1:
        ord += 1
        cycle_elem = sr ** ord % p
    
    return ord

def order_and_primitivity_of_t(N: int) -> dict:
    p = 23
    t = 12

    Zp = [i for i in range(1, p) if gcd(i, p) == 1]

    if t not in Zp:
        return None
    
    ord_t = -1

    for ord in range(1, len(Zp)):
        if t ** ord % p == 1:
            ord_t = ord
            break
    
    if is_prime(p) and ord_t == p - 1:
        is_generator = True
    else:
        is_generator = False

    return {
        "order": ord_t,
        "is_generator": is_generator
    }

def divs(n):
    divs = []
    for i in range(1, int(n ** 0.5) + 1):
        if n % i == 0:
            divs.append(i)
            divs.append(n // i)
    return divs

def generators_of_Zm_star(N: int) -> list:
    """
    Порядок пораждающего элемента равен порядку группы, потому он в степени порядка группы должен давать нейтральный элемент, 
    при этом он же в степени какого-либо своего делителя не должен давать нейтральный
    """
    m = 3 + N % 20

    Zm = [i for i in range(1, m) if gcd(i, m) == 1]
    generators = []

    for i in Zm:
        if i ** len(Zm) % m == 1:
            div = divs(i) 
            for j in div:
                if i ** j % m == 1:
                    break
            else:
                generators.append(i)

    return generators

def cyclic_subgroup_in_Zm_additive(N: int) -> dict:
    """
    Умножаем элемент на i > 2 до тех пор, пока не получим этот же элемент
    Элемент циклической группы является пораждающим, если НОД порядка группы и степени пораждающего элемента равен нулю
    """
    m = 4 + N % 20
    t = 17

    Zm = Zmod(m)

    el = t % m

    if el == 0:
        return {
            "Подгруппа": [el],
            "Пораждающие": [el]
        }

    subgroup = [el]

    generators = []
    el_div = {el: 1}

    i = 2
    while el * i % m != el:
        subgroup.append(el * i % m)
        el_div[el * i % m] = i
        i += 1
    
    for i in subgroup:
        if gcd(len(subgroup), el_div[i]) == 1:
            generators.append(i)

    return {
        "Подгруппа": sorted(subgroup),
        "Пораждающие": sorted(generators)
    }

def isomorphism_of_cyclic_subgroup_Zm_star(N: int) -> dict:
    """
    Любая конечная циклическая группа изоморфна аддитивной группе остатков, она же в свою очередь изоморфна циклической подгруппе 
    Sn, порожденной (1234...n)
    """
    m = 3 + N % 20
    t = 12
    Zm = Integers(m).unit_group()
    el = t % m

    subgroup = []

    i = 2
    while el ** i % m != el:
        subgroup.append(el ** i % m)
        i += 1

    subgroup.sort()
    
    return {
        "Подгруппа": subgroup,
        "Изоморфна": "<(" + "".join(f"{i}" for i in range(1, len(subgroup) + 1)) + ")>"
    }

def polinoms(N) -> dict:
    answer1 = []
    answer2 = []

    for a in range(4):
        f = a ** 9
        for i in range(9):
            f += ((i + N) % 4) * a ** i
        
        if f % 4 == 0:
            answer1.append(a)
    
    for b in range(7):
        f = 0
        for i in range(7):
            f += ((i + N) % 7) * b ** i
        
        if f % 7 == 0:
            answer2.append(b)

    return {
        "polinom 1": answer1,
        "polinom 2": answer2
    }

def factorization_polinoms(N) -> dict:
    c = [(i + N) % 5 for i in range(5)]
    d = [(i + N) % 9 for i in range(9)]

    multiples2 = []

    F5 = GF(5)
    R5 = PolynomialRing(F5, "x")
    x = R5.gen()

    f5 = x ** 5 + sum(c[i] * x ** i for i in range(5))
    multiples1 = f5

    if not f5.is_irreducible():
        multiples1 = f5.factor()
    
    F9 = GF(9)
    R9 = PolynomialRing(F9, "x")
    x = R5.gen()

    f9 = x ** 4 + sum(d[i] * x ** i for i in range(4))
    multiples2 = f9

    if not f9.is_irreducible():
        multiples2 = f9.factor()
    
    return {
        "Разложение 1": multiples1,
        "Разложение 2": multiples2
    }

def gcd_polinoms(N) -> dict:
    r = [(i + N) % 11 for i in range(11)]
    s = [(i + N) % 11 for i in range(11)]

    F11 = GF(11)
    R = PolynomialRing(F11, "x")
    x = R.gen()

    f = sum(r[i] * x ** i for i in range(8))
    g = sum(s[i] * x ** i for i in range(4))

    gcd_pol, u, v = xgcd(f, g)

    return {
        "Нод": gcd_pol,
        "Коэффиценты Безу": (u, v)
    }

def reverse_polinom(N) -> str:
    """
    Пусть есть некоторый полином u(x), что u*f = 1 (mod g). Тогда для любого v: u*f + v*g = 1 (mod g). Это есть тождество Безу.
    """

    s = [(i + N) % 11 for i in range(11)]

    F13 = GF(13)
    R = PolynomialRing(F13, "x")
    x = R.gen()

    f = s[2] * x ** 2 + s[1] * x + s[0]
    g = x ** 8 + x ** 4 + x ** 3 + 6 * x + 2

    if not gcd(f, g) == 1:
        return "Нет обратного"
    
    gcd_pol, u, v = xgcd(f, g)

    return str(u)

def irreducible_polynoms(q: int, d: int) -> list:
    """
    Находим все перестановки элементов поля, строим полиномы, проверяем на приводимость
    """
    if d < 1:
        return []
    
    F = GF(q)
    R = PolynomialRing(F, 'x')
    x = R.gen()

    irreducibles = []

    for coeffs in product(F, repeat=d + 1):
        if coeffs[-1] == 0:
            continue
        f = sum(coeffs[i] * x**i for i in range(d + 1))
        if f.is_irreducible():
            irreducibles.append(f)
    return irreducibles
