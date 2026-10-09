"""Exception hierarchy for bolzakit."""


class BolzakitError(Exception):
    """Base class for all bolzakit errors."""


class OutOfValidityRange(BolzakitError, ValueError):
    """Input lies outside the region where the method is validated."""


class DataUnavailable(BolzakitError, RuntimeError):
    """Required data could not be obtained or failed its checksum."""


class OptionalDependencyMissing(BolzakitError, ImportError):
    """An optional dependency (pip extra) is not installed."""


class CertificationFailed(BolzakitError):
    """A certified computation did not separate its enclosures (raised only when strict=True)."""
