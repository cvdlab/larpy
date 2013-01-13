# algoritmo estrazione faccette
# =============================

# input:  matrice caratteristica $M_d$;
# output:  matrice caratteristica $M_{d-1}$.

# decomposizione iniziale
# -----------------------

# 1. Calcolo faccette come sottoinsiemi $S_{i,j}$ di vertici comuni alle coppie di celle $S_i, S_j \in \Lambda_d(X)$ ($i,j\in \N$). $\sharp S_{i,j} \geq \dim(S_i) = \dim(S_j) = d$


# 2. Partizione delle faccette di cardinalità $> d-1$ in sottoinsiemi "piatti", ovvero con dimensione del guscio affine pari a $d-1$.


# 3. Sia $M_{d-1}$ la matrice caratteristica dei sottoinsiemi $S_{i,j}$, ordinata per righe rispetto al numero (crescente) di vertici incidenti. Si vuole sostituire ad ogni riga con somma maggiore di $d$ un insieme di righe di somma $d$ corrispondenti ad un insieme connesso di faccette piatte con gli stessi vertici, ovvero con la stessa somma sulle colonne.


# 4.  Questa decomposizione delle righe è unica? probabilmente no. Dovendo scegliere tra più decomposizioni ammissibili, sceglieremmo quella (unica?) dove la catena cobordo del bordo abbia---a tratti---lo stesso guscio affine delle faccette esterne adiacenti al bordo.


# decomposizione fine
# -------------------

# 5. Consideriamo la prima faccetta di cardinalità maggiore di $d$, sia $S_k$, di indice $k$ ed estraiamo la sottomatrice massimale $N_k$ delle righe di $M_{d-1}$ con indice di riga minore di $k$  e somma per colonne pari a $2M_k$, ovvero delle faccette adiacenti di cardinalità eguale a $d$.

# 6. Partizioniamo $M_k$ (utilizzando $N_k$ e/o $N_k^T.N_k$ ed $N_k.N_k^T$) in più righe $R_k$, che sostituiamo ad $M_k$.

# 7. Incrementiamo $k$ e ritorniamo al punto 5, fino a terminare.
