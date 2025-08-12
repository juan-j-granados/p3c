"""
p3c - Cubic polynomial transformation based on the complex plane
for regularizing nearly singular integrals in the boundary element method.
"""

__author__ = "Juan J. Granados, Rafael Gallego"
__version__ = "1.0.0"
__license__ = "MIT"

from .core import build_p3c

__all__ = ["build_p3c"]