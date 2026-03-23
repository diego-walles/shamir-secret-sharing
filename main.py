from typing import List, Tuple


def mod_inverse(a: int, p: int) -> int:
    """
    Retorna el inverso multiplicativo de a modulo p.
    Requiere que gcd(a, p) = 1.
    """
    a %= p
    if a == 0:
        raise ValueError("No existe inverso modular para 0.")
    return pow(a, -1, p)


def build_vandermonde(points: List[Tuple[int, int]], p: int) -> Tuple[List[List[int]], List[int]]:
    """
    Construye el sistema lineal modular:
        A * c = b (mod p)
    donde:
        - A es la matriz de Vandermonde
        - c = [a0, a1, ..., a_{n-1}]
        - b = [y0, y1, ..., y_{n-1}]
    Si hay n puntos, se asume un polinomio de grado n-1.
    """
    n = len(points)
    A = []
    b = []

    for x, y in points:
        row = []
        value = 1
        for _ in range(n):
            row.append(value % p)
            value = (value * x) % p
        A.append(row)
        b.append(y % p)

    return A, b


def solve_modular_system(A: List[List[int]], b: List[int], p: int) -> List[int]:
    """
    Resuelve el sistema A*x = b (mod p) usando eliminación de Gauss-Jordan.
    """
    n = len(A)

    # Matriz aumentada
    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    for col in range(n):
        # Buscar pivote no nulo
        pivot = None
        for row in range(col, n):
            if M[row][col] % p != 0:
                pivot = row
                break

        if pivot is None:
            raise ValueError("La matriz no es invertible modulo p.")

        # Intercambiar filas si es necesario
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]

        # Normalizar fila pivote
        inv_pivot = mod_inverse(M[col][col], p)
        for j in range(col, n + 1):
            M[col][j] = (M[col][j] * inv_pivot) % p

        # Eliminar en las otras filas
        for row in range(n):
            if row != col and M[row][col] % p != 0:
                factor = M[row][col] % p
                for j in range(col, n + 1):
                    M[row][j] = (M[row][j] - factor * M[col][j]) % p

    # Extraer solución
    solution = [M[i][n] % p for i in range(n)]
    return solution


def polynomial_to_string(coeffs: List[int], p: int) -> str:
    """
    Convierte la lista de coeficientes [a0, a1, a2, ...]
    en una representación legible del polinomio modulo p.
    """
    terms = []
    for i, c in enumerate(coeffs):
        c %= p
        if i == 0:
            terms.append(f"{c}")
        elif i == 1:
            terms.append(f"{c}*x")
        else:
            terms.append(f"{c}*x^{i}")
    return " + ".join(terms) + f"  (mod {p})"


def evaluate_polynomial(coeffs: List[int], x: int, p: int) -> int:
    """
    Evalúa el polinomio en x modulo p.
    coeffs = [a0, a1, a2, ...]
    """
    result = 0
    power = 1
    for c in coeffs:
        result = (result + c * power) % p
        power = (power * x) % p
    return result


def reconstruct_secret(points: List[Tuple[int, int]], p: int) -> Tuple[List[int], int]:
    """
    Reconstruye el polinomio usando todos los puntos dados
    y retorna (coeficientes, secreto).
    El secreto es a0 = q(0).
    """
    if len(points) == 0:
        raise ValueError("Debe proporcionar al menos un punto.")

    xs = [x for x, _ in points]
    if len(xs) != len(set(xs)):
        raise ValueError("Los valores x deben ser distintos.")

    A, b = build_vandermonde(points, p)
    coeffs = solve_modular_system(A, b, p)
    secret = coeffs[0] % p
    return coeffs, secret


def main() -> None:
    # Ejemplo del ejercicio para poder despejarlo y para claridad, esto me debería salir en la terminal cuando yo ejecute el main.py 
    p = 1301
    points = [(2, 600), (5, 960), (7, 120), (1, 360)]

    print("=== Shamir Secret Sharing - Reconstrucción por matrices modulares ===")
    print(f"Primo p = {p}")
    print(f"Puntos dados = {points}")
    print()

    coeffs, secret = reconstruct_secret(points, p)

    print("Coeficientes del polinomio q(x):")
    for i, c in enumerate(coeffs):
        print(f"a{i} = {c}")

    print()
    print("Polinomio reconstruido:")
    print(f"q(x) = {polynomial_to_string(coeffs, p)}")

    print()
    print(f"Secreto recuperado S = q(0) = {secret}")

    print()
    print("Verificación con los puntos dados:")
    for x, y in points:
        calculated = evaluate_polynomial(coeffs, x, p)
        print(f"q({x}) mod {p} = {calculated} | esperado = {y} | OK = {calculated == (y % p)}")


if __name__ == "__main__":
    main()