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

// Exercicio 4 - 8.7
// Dados
Dbl M5 = 500000 // Massa [kg]
Dbl F5 = 1/M5 // Inverso da massa
Dbl K5 = 300000*10^3 // Rigidez [N/m]
Dbl Xi5 = 0.02 // Taxa de Amortecimento [%]
Dbl dt5 = 0.05 // Intervalo de tempo da funcao carregamento
Dbl N5 = 128 // Numero de passos
Dbl W5 = Sqrt[K5/M5] // Frequencia Natural [rad/s]
Dbl gama5 = dt5*N5 // Tempo de observacao [s]

//matriz do carregamento:

Dbl P5 =[0.00 0.00
	 0.05 1545084.97
	 0.10 2938925.26
	 0.15 4045084.97
	 0.20 4755282.58
	 0.25 5000000.00
	 0.30 4755282.58
	 0.35 4045084.97
	 0.40 2938925.26
	 0.45 1545084.97
	 0.50 0
	 gama5 0]

//Resolucao:

Dbl qF = Resample(P5,N5,dt5) // Reordena q para cada intervalo de tempo dt
Dbl FFTP = FFT_Real(qF,dt5) // Retorna a Transformada Rapida de Fourier de um sinal real
Dbl FH = FSolve(FFTP,W5,Xi5,F5,dt5) // Retorna o produto da funcao da forca no dominio de frequencia pela funcao de resposta no dominio da frequencia
Dbl U5 = iFFT_Real(FH,dt5) // Retorna a Inversa da Transformada Rapida de Fourier de um sinal complexo
	