# v0.2 Detailed Derivations

## A. The second derivative coupling identity

Define the derivative-coupling matrix

$$
d_{ab}=\langle\phi_a|\partial\phi_b\rangle.
$$

Differentiate:

$$
\partial d_{ab}
=
\langle\partial\phi_a|\partial\phi_b\rangle
+
\langle\phi_a|\partial^2\phi_b\rangle.
$$

Insert completeness,

$$
I=\sum_c|\phi_c\rangle\langle\phi_c|,
$$

into the first term:

$$
\langle\partial\phi_a|\partial\phi_b\rangle
=
\sum_c
\langle\partial\phi_a|\phi_c\rangle
\langle\phi_c|\partial\phi_b\rangle.
$$

From differentiated orthonormality,

$$
\langle\partial\phi_a|\phi_c\rangle=-d_{ac}.
$$

Therefore

$$
\langle\partial\phi_a|\partial\phi_b\rangle
=
-\sum_c d_{ac}d_{cb}.
$$

Hence

$$
\boxed{
\tau_{ab}
=
\langle\phi_a|\partial^2\phi_b\rangle
=
\partial d_{ab}
+
\sum_c d_{ac}d_{cb}.
}
$$

or

$$
\boxed{\tau=d'+d^2.}
$$

## B. Expansion of the covariant kinetic operator

Apply

$$
(\partial I+d)^2\chi.
$$

First,

$$
(\partial I+d)\chi=\chi'+d\chi.
$$

Apply it again:

$$
\partial(\chi'+d\chi)+d(\chi'+d\chi)
$$

$$
=
\chi''+d'\chi+d\chi'+d\chi'+d^2\chi.
$$

Therefore,

$$
\boxed{
(\partial I+d)^2\chi
=
\chi''
+
2d\chi'
+
(d'+d^2)\chi.
}
$$

This reproduces the Born-Huang kinetic coupling exactly within a complete electronic
subspace.

## C. Diabatic-to-adiabatic operator transformation

Let

$$
\psi_d=U\chi.
$$

Differentiate:

$$
\partial(U\chi)=U'\chi+U\chi'.
$$

Because

$$
d=U^\dagger U',
$$

we have

$$
U'=Ud.
$$

Thus

$$
\partial(U\chi)
=
U(d\chi+\chi')
=
U(\partial+d)\chi.
$$

Differentiate again:

$$
\partial^2(U\chi)
=
U(\partial+d)^2\chi.
$$

Therefore,

$$
U^\dagger
\left(
-\frac{1}{2M}\partial^2 I+V_d
\right)
U
=
-\frac{1}{2M}(\partial+d)^2+E.
$$

This is the exact representation equivalence.

## D. Energy-conserving spawn momentum

Require

$$
K_a+E_a=K_b+E_b.
$$

In 1D,

$$
\frac{p_a^2}{2M}+E_a
=
\frac{p_b^2}{2M}+E_b.
$$

Multiply by $2M$:

$$
p_a^2+2ME_a=p_b^2+2ME_b.
$$

Hence

$$
\boxed{
p_b^2=p_a^2+2M(E_a-E_b).
}
$$

The existence of a real child under this local rule requires

$$
p_a^2+2M(E_a-E_b)\ge0.
$$
