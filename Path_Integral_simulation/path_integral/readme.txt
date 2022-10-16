Test per il codice harmonic_oscillator_simple.f
che implementa l'oscillatore armonico semplice

parametri: N = 10, eta = 0.1, delta_metropolis = 0.5
           10 milioni di misure, una presa ogni 10 spazzate del reticolo
           scarto il primo milione di misure
           <x^2> = 1.083(3)
           tempo di autocorrelazione: fra 25 e 26 (in unita` di distanza
           fra una misura e l'altra).
 
===========================================================================


Nel codice dboson.f abbiamo implementato una versione semplificata
del codice per i due bosoni in potenziale armonico in cui, alla fine
di ogni spazzata di metropolis locale, viene tentato un cambio
delle condizioni al contorno, cioe` viene implementato
l'algoritmo visto a lezione scegliendo sempre j_0 = N.
La cosa non e` ottimale ma funziona ragionevolmente bene.

I due vettori corrispondenti ai due cammini sono implementati
con un unico vettore lungo 2N, le condizioni al contorno 
sono implementate nella funzione "geometry" che viene 
aggiornata alla fine di ogni metropolis per implementare 
l'eventuale cambio di condizioni al contorno.

  
Riportiamo il seguente benchmark:

N = 30, eta = 0.1, delta = 0.5
10^6 misure, ognuna presa ogni 100 spazzate, 10^5 spazzate iniziali eliminate per termalizzazione:

se indichiamo con x ed y i due cammini:

<(x^2 + y^2)/2>  = 0.5289(7)    (valore atteso nel limite continuo  0.52868276)
 
< segno > = 0.0497(14)          (valore atteso nel limite continuo  0.04978707)



