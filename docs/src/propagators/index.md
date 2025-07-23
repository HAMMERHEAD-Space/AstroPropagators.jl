# Propagation Methods

AstroPropagators.jl provides a comprehensive suite of orbital propagation methods, each optimized for different orbital regimes and applications. This section provides detailed documentation for each propagation method.

## Available Methods

### [Cowell Method](cowell.md)

The Cowell method propagates orbits using Cartesian coordinates and is the most widely used approach for orbital mechanics. It provides excellent accuracy and can handle any type of force model.

### [Gauss Variational Equations](gauss_ve.md)

Gauss Variational Equations propagate orbital elements directly, making them ideal for studying how perturbations affect orbital parameters over time.

### [EDromo Method](edromo.md)

EDromo uses regularized elements and fictitious time to handle challenging orbital regimes where traditional methods struggle.

### [Kustaanheimo-Stiefel](kustaanheimo_stiefel.md)

The K-S method uses 4D regularized coordinates to eliminate singularities and provide numerical stability for extreme orbital conditions.

### [Milankovich Elements](milankovich.md)

Milankovich elements use angular momentum and eccentricity vectors to provide a non-singular alternative to classical elements.

### [Stiefel-Scheifele](stiefel_scheifele.md)

The Stiefel-Scheifele method provides canonical regularization for orbital dynamics with excellent numerical properties.

### [Unified State Model](usm.md)

USM provides a unified framework with three variants (USM7, USM6, USMEM) that combine velocity hodograph components with different attitude parameterizations.