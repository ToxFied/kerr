# Physics

Units are geometric, $G = c = M = 1$, so lengths are in units of the hole's mass $M$. The spin is $a = J/M$, with $|a| \le 0.998$. Everything below is implemented in [`src/physics.h`](../src/physics.h).

## Spacetime

The Kerr metric in Boyer-Lindquist coordinates $(t, r, \theta, \phi)$ uses

$$
\Sigma = r^2 + a^2\cos^2\theta, \qquad
\Delta = r^2 - 2r + a^2, \qquad
A = (r^2 + a^2)^2 - a^2 \Delta \sin^2\theta .
$$

The outer horizon is at $r_+ = 1 + \sqrt{1 - a^2}$. Circular photon orbits in the equatorial plane sit at

$$
r_{\rm ph} = 2\left[1 + \cos\left(\tfrac{2}{3}\arccos(\mp a)\right)\right],
$$

the upper sign for orbits that co-rotate with the hole. Between the two lies the photon shell, where spherical photon orbits exist.

## Null geodesics

A photon has three constants of motion besides its rest mass of zero: energy $E$, angular momentum about the spin axis $L$, and Carter's constant $Q$. Only the ratios $\ell = L/E$ and $\eta = Q/E^2$ matter for the path.

Carter's constant separates the equations. In Mino time $\lambda$, defined by $d\lambda = d\sigma/\Sigma$ with $\sigma$ the affine parameter, the radial and polar motions decouple:

$$
\left(\frac{dr}{d\lambda}\right)^2 = R(r) = (r^2 + a^2 - a\ell)^2 - \Delta\left[\eta + (\ell - a)^2\right],
$$

$$
\left(\frac{d\mu}{d\lambda}\right)^2 = M(\mu) = \eta - (\eta + \ell^2 - a^2)\,\mu^2 - a^2\mu^4, \qquad \mu = \cos\theta,
$$

$$
\frac{d\phi}{d\lambda} = \frac{a\,(r^2 + a^2 - a\ell)}{\Delta} - a + \frac{\ell}{1 - \mu^2}.
$$

The first-order forms need the sign of each square root tracked through turning points. Differentiating once gives second-order equations that pass turning points on their own:

$$
\frac{d^2 r}{d\lambda^2} = \tfrac12 R'(r), \qquad \frac{d^2\mu}{d\lambda^2} = \tfrac12 M'(\mu).
$$

### Why the radial equation uses $u = 1/r$

Far from the hole $R(r) \approx r^4$, while the part that carries the impact parameter is $\sim b^2 r^2$, smaller by a factor $b^2/r^2$. At $r = 10^4$ and $b = 5$ that is $10^{-7}$, below what RK4 round-off preserves over a few hundred steps, and rays turned around at the wrong radius. Writing $u = 1/r$ and $U(u) = u^4 R(1/u)$ fixes the conditioning:

$$
U(u) = (1 + c\,u^2)^2 - K u^2 (1 - 2u + a^2 u^2), \qquad c = a^2 - a\ell, \quad K = \eta + (\ell - a)^2,
$$

$$
\frac{d^2u}{d\lambda^2} = \tfrac12 U'(u) = 2cu\,(1 + cu^2) - Ku\,(1 - 3u + 2a^2u^2).
$$

Now the impact-parameter term is order one at every radius. At $a = 0$ this is Binet's equation for light, $u'' + u = 3u^2$, in a rescaled parameter.

### Integration

The state $(u, \dot u, \mu, \dot\mu, \phi)$ is advanced with fourth-order Runge-Kutta. The step is chosen so that each step changes $r$ by about $h\,r$ and $\mu$ and $\phi$ by about $h$, with $h = 0.02$ for stills. The viewer's preview passes let $h$ grow with $r$ beyond $20\,M$, where rays are nearly straight.

A ray hits the disk when $\mu$ changes sign between steps at a radius inside the disk; the crossing is found by bisecting the RK4 sub-step to $2^{-24}$ of its length. It falls into the hole once it reaches $r < \tfrac12 (r_+ + r_{\rm ph,min})$: turning points of rays arriving from outside lie beyond the photon shell, so a ray inside that radius cannot come back. It escapes to the sky once it passes $r = 10^4$ moving outward, and its direction there comes from the Mino-time rates converted to Cartesian velocity.

Each crossing of the equatorial plane outside the disk raises the ray's image order: order 0 is the direct image of the disk, 1 the secondary image, and 2 or more the photon ring.

## Camera

The camera is a zero-angular-momentum observer (ZAMO), whose local frame is

$$
e^{(t)} = \sqrt{\tfrac{\Sigma\Delta}{A}}\,dt, \quad
e^{(r)} = \sqrt{\tfrac{\Sigma}{\Delta}}\,dr, \quad
e^{(\theta)} = \sqrt{\Sigma}\,d\theta, \quad
e^{(\phi)} = \sqrt{\tfrac{A}{\Sigma}}\sin\theta\,(d\phi - \omega\,dt),
$$

with frame-dragging rate $\omega = 2ar/A$. A pixel defines a viewing direction $\mathbf n$ in this frame. The photon arriving there travels along $\mathbf d = -\mathbf n$, and with locally measured energy 1 its momentum components are

$$
p_\phi = L = \sqrt{\tfrac{A}{\Sigma}}\sin\theta\, d_\phi, \qquad
E = \sqrt{\tfrac{\Sigma\Delta}{A}} + \omega L, \qquad
p_\theta = \sqrt{\Sigma}\, d_\theta, \qquad
p_r = \sqrt{\tfrac{\Sigma}{\Delta}}\, d_r .
$$

Dividing by $E$ gives $\ell$, and Carter's constant follows as $\eta = p_\theta^2 + \cos^2\theta\,(\ell^2/\sin^2\theta - a^2)$. The initial Mino-time rates are $\dot r = \Delta p_r$ and $\dot\mu = -\sin\theta\, p_\theta$. Both are negated because the ray is traced backward in time.

## Accretion disk

The disk is thin, optically thick, and lies in the equatorial plane between the innermost stable circular orbit (ISCO) and $r_{\rm out}$. Its gas follows circular geodesics (Bardeen, Press & Teukolsky 1972):

$$
\Omega = \frac{1}{r^{3/2} + a}, \qquad
u^t = \frac{r^{3/2} + a}{r^{3/4}\sqrt{r^{3/2} - 3r^{1/2} + 2a}},
$$

$$
E = \frac{r^{3/2} - 2r^{1/2} + a}{r^{3/4}\sqrt{r^{3/2} - 3r^{1/2} + 2a}}, \qquad
L = \frac{r^2 - 2ar^{1/2} + a^2}{r^{3/4}\sqrt{r^{3/2} - 3r^{1/2} + 2a}} .
$$

A negative $a$ describes a disk that orbits against the spin. The ISCO is

$$
r_{\rm isco} = 3 + Z_2 \mp \sqrt{(3 - Z_1)(3 + Z_1 + 2Z_2)},
$$

$$
Z_1 = 1 + (1 - a^2)^{1/3}\left[(1 + a)^{1/3} + (1 - a)^{1/3}\right], \qquad Z_2 = \sqrt{3a^2 + Z_1^2}.
$$

The radiative efficiency, the fraction of rest mass released by gas spiralling down to the ISCO, is $1 - E(r_{\rm isco})$: 5.7% at $a = 0$, 15.6% at 0.9 and 32.1% at 0.998.

The flux leaving each face of the disk is the Page-Thorne (1974) result for a stationary, thin disk with no torque at the ISCO:

$$
F(r) = -\frac{\dot M}{4\pi r}\,\frac{\Omega_{,r}}{(E - \Omega L)^2}\int_{r_{\rm isco}}^{r} (E - \Omega L)\,L_{,r}\,dr ,
$$

evaluated numerically for the current spin. The disk radiates as a local blackbody with $T \propto F^{1/4}$, scaled so its peak equals the chosen temperature.

## Redshift and color

For a photon with $\ell$ arriving at the camera, light from disk radius $r$ is shifted by

$$
g = \frac{\nu_{\rm obs}}{\nu_{\rm emit}} = \frac{1}{E\,u^t\,(1 - \Omega\,\ell)},
$$

where $E$ is the photon energy at infinity when the camera measures energy 1. This one factor includes gravitational redshift, the Doppler shift of the orbiting gas, and the camera's own frame.

Since $I_\nu / \nu^3$ is conserved along a ray, a Planck spectrum at temperature $T$ arrives as a Planck spectrum at $gT$. The observed spectrum is integrated against the CIE 1931 color matching functions (multi-lobe fit of Wyman, Sloan & Shirley 2013) to get XYZ, then converted to linear sRGB with a D65 white point. A 6504 K blackbody comes out white.

## Sky

Stars come from the Yale Bright Star Catalog (Hoffleit & Warren 1991): 9,096 stars to visual magnitude 6.5. Each star's temperature comes from its B-V color (Ballesteros 2012),

$$
T = 4600\,\mathrm{K}\left(\frac{1}{0.92\,(B-V) + 1.7} + \frac{1}{0.92\,(B-V) + 0.62}\right),
$$

and its brightness from its magnitude, $10^{-0.4V}$. Starlight falling to the camera is blueshifted by $1/\alpha$, where $\alpha = \sqrt{\Sigma\Delta/A}$ is the observer's lapse. The frame-dragging correction $\omega L$ is below 0.1% at the default distance and is left out. Each star is drawn as a narrow Gaussian in its true sky direction plus a faint halo three times wider, like a camera's point-spread function. Because the shapes live on the sky, lensing stretches magnified images into arcs and squeezes demagnified ones.

The catalog is rotated so that Orion (RA 5h35m, Dec -1.5°) lies behind the hole for the default camera, with celestial north up and east to the left, as seen from Earth.

## Shadow edge

The critical curve drawn in the viewer is computed for the actual camera rather than from a far-observer formula. For 180 directions around the view center, bisection on the screen radius separates rays that are captured from rays that escape (with the disk removed). The bisection takes about 10 ms per view.

## Display

The disk's brightness spans several orders of magnitude, especially at high spin, where the Doppler-boosted inner edge outshines the outer disk. The display uses an asinh stretch (Lupton et al. 2004),

$$
L' = \frac{\operatorname{asinh}(\beta L)}{\operatorname{asinh}\beta}, \qquad \beta = 8,
$$

applied to luminance, with all three channels scaled by $L'/L$ so colors keep their hue. Exposure is fixed: disk luminance is normalized so the peak temperature at $g = 1$ has $Y = 1$, so stars and disk keep their real relative brightness from any viewpoint.

## Checks

`./blackhole --check` runs these comparisons:

- At $a = 0$, rays from infinity are captured exactly when $b < 3\sqrt3$.
- At $a = 0$ and $b = 100$, the bending angle matches $\alpha = \frac{4}{b} + \frac{15\pi}{4b^2} + \frac{128}{3b^3} + \frac{3465\pi}{64b^4} + \frac{3584}{5b^5}$ (Keeton & Petters 2005).
- At $a = 0.5$ and $0.9$, an equatorial observer sees the shadow's left and right edges at the $\ell$ of the co- and counter-rotating circular photon orbits, $\ell = -(r^3 - 3r^2 + a^2 r + a^2)/[a(r - 1)]$ (Bardeen 1973).
- The closed-form ISCO agrees with the minimum of $E(r)$.
- A ray that circles the photon shell keeps $\dot u^2 = U(u)$.
- From $10^5\,M$, both images of a star two Einstein radii off-axis solve $\beta = \theta - \alpha(b)$ with $b = D\sin\theta$.
- At $a = 0$ the disk efficiency is $1 - \sqrt{8/9}$ and its flux peaks near $9.55\,M$.
- Betelgeuse appears up and to the left of Rigel, and the catalog's aim point lies exactly behind the hole.
