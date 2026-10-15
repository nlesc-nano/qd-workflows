# IR and Raman

Intensities from the CP2K finite-field single points (`field` step of the chain), peak positions from the MACE-MH-1 Hessian; assembled by `vibspec` in the CPU stage.

## From the finite-field single points to spectra

Thirteen single points in a uniform field $\mathbf E$ (amplitude $E = 0.002$ a.u.), restarted from the
geo_opt wavefunction. The forces respond as

$$
Z^*_{i,ab} = \frac{\partial F_{ib}}{\partial E_a} \approx \frac{F_{ib}(+E\hat a) - F_{ib}(-E\hat a)}{2E}, \qquad
\frac{\partial \alpha_{ac}}{\partial R_{ib}} = \frac{\partial^2 F_{ib}}{\partial E_a\,\partial E_c}
$$

$$
\frac{\partial^2 F}{\partial E_a \partial E_c} \approx \frac{F(E,E) + F(-E,-E) - F(E,0) - F(-E,0) - F(0,E) - F(0,-E) + 2F(0,0)}{2E^2}
$$

Projected on the MACE normal modes $L_k$ (mass-weighted, masses $m_i$):

$$
\frac{\partial \mu_a}{\partial Q_k} = \sum_{i,b} Z^*_{i,ab}\,\frac{L_{ib,k}}{\sqrt{m_i}}, \qquad
\frac{\partial \alpha_{ac}}{\partial Q_k} = \sum_{i,b} \frac{\partial \alpha_{ac}}{\partial R_{ib}}\,\frac{L_{ib,k}}{\sqrt{m_i}}
$$

$$
A_k = \frac{N_A \pi}{3c}\left|\frac{\partial \boldsymbol\mu}{\partial Q_k}\right|^2, \qquad
S_k = 45\,a'^2 + 7\,\gamma'^2, \qquad \rho_k = \frac{3\gamma'^2}{45 a'^2 + 4\gamma'^2}
$$

$$
a' = \tfrac13\,\mathrm{tr}\,\alpha', \qquad
\gamma'^2 = \tfrac12\left[(\alpha'_{xx} - \alpha'_{yy})^2 + (\alpha'_{yy} - \alpha'_{zz})^2 + (\alpha'_{zz} - \alpha'_{xx})^2
+ 6(\alpha'^2_{xy} + \alpha'^2_{yz} + \alpha'^2_{zx})\right]
$$

Peak positions are the harmonic MACE-MH-1 frequencies.

