// Dinâmica das Estruturas - Lista 4
// Exercicio 2 - 7.1

// Dados
Dbl K1 = 2500000 // Rigidez total [N/m]
Dbl M1 = 17500 // Massa [kg]
Dbl W1 = Sqrt(K1/M1) // Frequência Natural [rad/s]
Dbl Xi1 = 0.00 // Amortecimento [%]

// Para dt = .01
Dbl dt1 = 0.0025
Dbl P01 = [0
68750
137500
206250
275000
343750
412500
481250
550000
481250
412500
343750
275000
206250
137500
68750
0]


// Calculando a integral de Duhamel
Dbl P1 = Zero(1000,1)
P1(1,1) = P01

Dbl Urec1 = Duhamel_Rectangle(P1, M1, W1, Xi1, dt1) // por integrao numerica de retangulos
Dbl Utrap1 = Duhamel_Trapeze(P1, M1, W1, Xi1, dt1) // de trapezeios
Dbl Usimp1 = Duhamel_Simpson(P1, M1, W1, Xi1, dt1) // por Simpson

// Exercicio 3 - 7.2
// Dados 
Dbl K2 = 7.5e6 // Rigidez [N/m]
Dbl M2 = 750e3 // Massa [kg]
Dbl W2 = Sqrt(K2/M2) // Frequência Natural [rad/s]
Dbl Xi2 = 0.0 // Amortecimento [%]

Dbl dt2 = 0.025
Dbl P02 = [0
80000
158000
223000
250000
125000
35000
7000
0]

Dbl P2 = Zero(800,1)
P2(1,1) = P02

Dbl Urec2 = Duhamel_Rectangle(P2, M2, W2, Xi2, dt2) // por integrao numerica de retangulos
Dbl Utrap2 = Duhamel_Trapeze(P2, M2, W2, Xi2, dt2) // de trapezeios
Dbl Usimp2 = Duhamel_Simpson(P2, M2, W2, Xi2, dt2) // por Simpson