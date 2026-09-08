---

# J. McLachlan TDVP derivation with real parameters

Let $$|\Psi\rangle=|\Psi(\theta_1,\ldots,\theta_M)\rangle$$ with real $\theta_\mu$.

Define $$|D_\mu\rangle=\partial_{\theta_\mu}|\Psi\rangle.$$

Then $$|\dot\Psi\rangle=\sum_\nu |D_\nu\rangle\dot\theta_\nu.$$

The TDSE residual is $$|r\rangle=i\sum_\nu|D_\nu\rangle\dot\theta_\nu-\hat H|\Psi\rangle.$$

Minimize $$\mathcal R=\langle r|r\rangle$$ with respect to each real $\dot\theta_\mu$.

Because $$\frac{\partial r}{\partial\dot\theta_\mu}=i|D_\mu\rangle,$$ stationarity gives $$0=2\operatorname{Re}\left\langle iD_\mu\middle|r\right\rangle.$$

Substitute $r$:

$$0=2\operatorname{Re}\left[\sum_\nu\langle iD_\mu|iD_\nu\rangle\dot\theta_\nu-\langle iD_\mu|H|\Psi\rangle\right].$$

Now $$\langle iD_\mu|iD_\nu\rangle=\langle D_\mu|D_\nu\rangle,$$ and $$\operatorname{Re}\langle iD_\mu|H|\Psi\rangle=\operatorname{Im}\langle D_\mu|H|\Psi\rangle.$$

Therefore,

$$\boxed{\sum_\nu\operatorname{Re}\langle D_\mu|D_\nu\rangle\dot\theta_\nu=\operatorname{Im}\langle D_\mu|H|\Psi\rangle.}$$

Or

$$\boxed{G\dot\theta=b.}$$

---
