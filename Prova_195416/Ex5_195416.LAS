// Carlos Henrique Chama Puga - 195416
// Prova - Dinamica das Estrturas - 2S2022
// Script LAS para construcao do espectro de resposta
// de uma ponte.

Input("ElCentro195416.LAS") // input do ElCentro

Dbl Dt = 0.02
Dbl Tmin = 0.00
Dbl Tmax = 3.00
Dbl ns = 150
Dbl Xi1 = 0.07
Dbl Tminps = 0

//Construcao do espectro verdadeiro  para estrutura com apoio:
Dbl {SD1, SV1, SA1} = Spectrum(Ag, Dt, Tmin, Tmax, ns, Xi1, "Lin")

// construcao do espectro verdadeiro para estrutura sem apoio
Dbl Xi2 = 0.03
Dbl {SD2, SV2, SA2} = Spectrum(Ag, Dt, Tmin, Tmax, ns, Xi2, "Lin")