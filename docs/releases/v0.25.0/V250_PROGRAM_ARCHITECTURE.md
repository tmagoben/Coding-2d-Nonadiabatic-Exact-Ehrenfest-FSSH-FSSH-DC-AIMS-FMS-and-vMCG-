# v0.25.0 program architecture

```text
              COMPLETE ARBITRARY-GEOMETRY OPERATOR PROVIDER
        q --> snapshot {H(q), K_a(q), D_a(q), M, state vectors}
         \                                      /
          \---- cross-snapshot overlap O_01 ---/
                              |
                 O_01 = U Sigma V^dagger
                              |
          +-------------------+-------------------+
          |                                       |
          v                                       v
 raw contraction evidence                 unitary polar W = U V^dagger
 singular retention                       coefficient frame transport
 condition number                                  |
 principal angle                                   |
          |                                         |
          +------------- fail-closed ---------------+
                              |
                              v
                RESTRICTED TDVP SYMMETRIC STEP

   F_n=-<c_n|K_n|c_n>             exp(-i H_n h/2)
           |                              |
           v                              v
     p_(n+1/2)                      W_01^dagger
           |                              |
           v                              v
 q_(n+1)=q_n+h M^-1 p_(n+1/2)      exp(-i H_(n+1) h/2)
           |                              |
           +------------+-----------------+
                        v
             F_(n+1)=-<c_(n+1)|K_(n+1)|c_(n+1)>
                        |
                        v
              p_(n+1)=p_(n+1/2)+h F_(n+1)/2
                        |
                        v
              fingerprint-bound step receipt
                        |
          +-------------+-------------------+
          |                                 |
          v                                 v
  validated v0.25.0                 deliberately closed
  analytic even/odd SOC             full multi-Gaussian TDVP
  complex-gauge covariance          adaptive widths / variable mass
  signed-step reversal              real PySCF SOC trajectories
  second-order convergence          ab-initio dynamics accuracy
```

## Receipt trust chain

Every step recomputes the stored SVD, polar factor, overlap metrics, force
expectations, Verlet endpoints, electronic Strang endpoint, energy, and time. The
trajectory then verifies step-to-step state continuity and freezes its scientific
claims in its serialized fingerprint.
