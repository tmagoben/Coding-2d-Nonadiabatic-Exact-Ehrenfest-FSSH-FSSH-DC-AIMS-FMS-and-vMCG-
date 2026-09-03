# v0.24.1 program architecture

```text
PINNED PYSCF 2.13.1                     SPIN ALGEBRA

real common-orbital SA-CASSCF           spin-pure roots (E,S,Mref)
  converged ROHF + CASSCF                 complete multiplet expansion
  CI roots + transition 1-RDMs            integer 2S/2M representation
  state-averaged spin-free 1-RDM           Clebsch-Gordan finite sum
              |                                      |
              |                     zero q=0 CG -----+
              |                        |             |
              |                   PySCF S+/S- ladder |
              +------------------------+-------------+
                                       v
PYSCF AO INTEGRALS              Wigner-reduced transition densities
  int1e_prinvxp                             |
  int2e_p1vxp1                              |
       |                                    |
       v                                    v
explicit BP-SOMF tensor ----------> direct state-interaction assembly
  J - 3/2 K_L - 3/2 K_R              H_sf, H_SOC, H_total
  one 1/(2c^2) prefactor             |root,S,M_S> order
  -i scalar-MO transform              time reversal + projectors
       |                                    |
       +------------------+-----------------+
                          v
                  STATIC SOC ADMISSION
                  - exact runtime/method identity
                  - direct-vs-JK SOMF cross-check
                  - Hermiticity + eigensystem
                  - complete multiplets
                  - time reversal + Kramers
                          |
            +-------------+------------------+
            |                                |
            v                                v
     fixed-geometry H_SOC              moving nuclei
        VALIDATED                     FAIL CLOSED
                                      no K_SOC
                                      no connection
                                      no cross-geometry overlap
```

The v0.24.0 OpenMolcas intake is inherited unchanged and remains an independent
cross-code route. Its native numeric cross-parser flag is not altered by PySCF static
SOC evidence.
