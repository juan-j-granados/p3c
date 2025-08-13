"""
p3c - Cubic polynomial transformation based on the complex plane
for regularizing nearly singular integrals in the boundary element method.

Authors:
    Juan J. Granados (ORCID: https://orcid.org/0000-0002-0707-2897)
    Rafael Gallego (ORCID: https://orcid.org/0000-0002-7260-0940)

License:
    MIT License (see LICENSE file for details)
"""

__author__ = "Juan J. Granados, Rafael Gallego"
__version__ = "1.0.0"
__license__ = "MIT"

from .core import build_p3c

__all__ = ["build_p3c"]