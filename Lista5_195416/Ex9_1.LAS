// Dynamic of Strucutres - IC956
// Segundo Semestre de 2022
// Exercicio 9.1 

Dbl T = 2
Dbl Dt1 = 0.05
Dbl Dt2 = 0.1
Dbl Dt3 = 0.2

Dbl P1 = Zero(T/Dt1+1,1)
Dbl P2 = Zero(T/Dt2+1,1)
Dbl P3 = Zero(T/Dt3+1,1)

Dbl M = 1
Dbl K = 2*Pi()*2*Pi()
Dbl W = Sqrt(K/M)
Dbl Xi = 0
Dbl C = 2*Xi*M*W

Dbl F0 = 0
Dbl U0 = 1
Dbl V0 = 0
Dbl A0 = -K*U0/M

Dbl {UNA1, VNA1, ANA1} = Newmark_Average_Acceleration(P1,M,C,K,U0,V0,A0,Dt1)
Dbl {UNA2, VNA2, ANA2} = Newmark_Average_Acceleration(P2,M,C,K,U0,V0,A0,Dt2)
Dbl {UNA3, VNA3, ANA3} = Newmark_Average_Acceleration(P3,M,C,K,U0,V0,A0,Dt3)

Dbl {UNL1, VNL1, ANL1} = Newmark_Linear_Acceleration(P1,M,C,K,U0,V0,A0,Dt1)
Dbl {UNL2, VNL2, ANL2} = Newmark_Linear_Acceleration(P2,M,C,K,U0,V0,A0,Dt2)
Dbl {UNL3, VNL3, ANL3} = Newmark_Linear_Acceleration(P3,M,C,K,U0,V0,A0,Dt3)