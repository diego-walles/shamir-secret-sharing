# shamir-secret-sharing algorithm - Modular Matrix Solution
Learn solution of nxn matrices, but with modular arithmetic Diego Walles - ANDES University COL

Implementación en Python del algoritmo de reconstrucción del secreto de Shamir en donde utilicé:

- aritmética modular
- resolución de sistemas lineales nxn
- eliminación de Gauss-Jordan módulo un primo `p`

## Descripción

Dado un conjunto de puntos `(x, y)` pertenecientes a un polinomio sobre `Z_p`, el programa:

1. construye la matriz de Vandermonde,
2. resuelve el sistema lineal módulo `p`,
3. recupera los coeficientes del polinomio,
4. obtiene el secreto como `q(0)`.

## Ejemplo usado conforme a indicaciones

Se probó con:

- `p = 1301`
- puntos:
  - `(2,600)`
  - `(5,960)`
  - `(7,120)`
  - `(1,360)`

## Resultado

El secreto recuperado es:

`S = 190`

Polinomio reconstruido:

`q(x) = 190 + 109x + 74x^2 + 1288x^3 (mod 1301)`

## Requisitos

- Python 3.8 o superior
- No utilicé librerías externas

## Ejecución

```bash
python main.py