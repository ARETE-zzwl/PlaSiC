# License and Attribution

Except for third-party materials stated otherwise, PlaSiC's C code, build and
test scripts, desktop application, this online tutorial, and the website source
are released under the **GNU General Public License version 3 or later**
(SPDX: `GPL-3.0-or-later`).

The GPL permits research, teaching, modification, redistribution, and
commercial use. When distributing binaries or derivative works governed by the
GPL, distributors must provide the complete corresponding source code, retain
applicable notices, and release GPL-covered derivative works under compatible
GPL terms. This software and documentation are provided without any warranty.

Third-party materials and external dependencies remain subject to their own
licenses and are not relicensed by being included in, or used with, this
project.

## Relationship to PlaSim

PlaSiC is a derivative of
[PlaSim (Planet Simulator)](https://github.com/HartmutBorth/PLASIM).
For scientific attribution, the upstream model description should be cited as
Fraedrich et al. (2005; full reference below).
Upstream PlaSim is licensed under `GPL-2.0-or-later`, and the original author's
copyright and license notices remain in force; a copy of the GPLv2 text is
included.
Since 2026, PlaSiC has reimplemented the model in C11, restructured the
software, removed the SimBA module, do the historical & 4xCO2 & 1pctCO2 experiments and added the desktop application and this
online tutorial. The combined PlaSiC work is released under `GPL-3.0-or-later`,
exercising the "GPL v2 or any later version" permission granted by the
upstream project.

## Acknowledgements

This project owes a great deal to the teaching materials and open-source work of
others. In particular, *Introduction to Climate Modelling* by Stocker (2011) and
the [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl) project
(Klöwer et al., 2024) guided much of the design of PlaSiC. The former shaped how
the physics and mathematics are presented, while the latter is an inspiring model
of an interactive, approachable atmospheric GCM. Without these two resources,
PlaSiC could not have been designed. The original authors are gratefully
acknowledged.

- Stocker, T. (2011). *Introduction to Climate Modelling*. Springer Science & Business Media.
- Klöwer, M., Gelbrecht, M., Hotta, D., Willmert, J., Silvestri, S., Wagner, G. L., White, A., Hatfield, S., Kimpson, T., Constantinou, N. C., & Hill, C. (2024). SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility. *Journal of Open Source Software*, *9*(98), 6323. [https://doi.org/10.21105/joss.06323](https://doi.org/10.21105/joss.06323)

