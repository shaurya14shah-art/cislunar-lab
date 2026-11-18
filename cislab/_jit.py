"""Optional Numba support.

Numba compiles the tight numerical loops (Monte Carlo propagation, N-body force
sums) to machine code.  If Numba is not installed everything still runs - just
slowly - because ``njit`` then becomes a do-nothing decorator.
"""
try:
    from numba import njit, prange  # type: ignore
    HAVE_NUMBA = True
except Exception:  # pragma: no cover - exercised only on machines without numba
    HAVE_NUMBA = False
    prange = range

    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]
        return lambda f: f
