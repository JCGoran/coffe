from abc import ABC, abstractmethod


class Representation(ABC):
    @abstractmethod
    def __init__(self, *args, **kwargs):
        pass

    def to_dict(self):
        """
        The representation of the class as a dictionary.
        """
        return {
            key: getattr(self, key)
            for key in dir(self.__class__)
            if hasattr(getattr(self.__class__, key), "__set__")
            and not key.startswith("__")
        }

    def __repr__(self):
        """
        User-friendly representation of the class
        """
        return f"<{self.__class__.__name__}({self.to_dict()})>"

    def __str__(self):
        """
        User-friendly representation of the class (string)
        """
        return f"{self.__class__.__name__}({self.to_dict()})"

    def _repr_html_(self):
        names = self.to_dict().keys()
        values = self.to_dict().values()
        temp = (
            ("<tr>" + ("<th>{}</th>" * len(names)).format(*names) + "</tr>")
            if names
            else ""
        )
        header = f"<thead>{temp}</thead>"

        body = "<tbody>" + ("<td>{}</td>" * len(values)).format(*values) + "</tbody>"

        return f"<table>{header}{body}</table>"


class Covariance(Representation):
    def __init__(
        self,
        *,
        r1: float,
        r2: float,
        l1: int,
        l2: int,
        z: float,
        value: float,
    ):
        self._r1 = r1
        self._r2 = r2
        self._l1 = l1
        self._l2 = l2
        self._z = z
        self._value = value

    @property
    def r1(self):
        return self._r1

    @property
    def r2(self):
        return self._r2

    @property
    def l1(self):
        return self._l1

    @property
    def l2(self):
        return self._l2

    @property
    def z(self):
        return self._z

    @property
    def value(self):
        return self._value


class AverageCovariance(Representation):
    def __init__(
        self,
        *,
        r1: float,
        r2: float,
        l1: int,
        l2: int,
        z_min: float,
        z_max: float,
        value: float,
    ):
        self._r1 = r1
        self._r2 = r2
        self._l1 = l1
        self._l2 = l2
        self._z_min = z_min
        self._z_max = z_max
        self._value = value

    @property
    def r1(self):
        return self._r1

    @property
    def r2(self):
        return self._r2

    @property
    def l1(self):
        return self._l1

    @property
    def l2(self):
        return self._l2

    @property
    def z_min(self):
        return self._z_min

    @property
    def z_max(self):
        return self._z_max

    @property
    def value(self):
        return self._value


class Corrfunc(Representation):
    def __init__(
        self,
        *,
        r: float,
        mu: float,
        z: float,
        value: float,
    ):
        self._r = r
        self._mu = mu
        self._z = z
        self._value = value

    @property
    def mu(self):
        return self._mu

    @property
    def r(self):
        return self._r

    @property
    def z(self):
        return self._z

    @property
    def value(self):
        return self._value


class Multipoles(Representation):
    def __init__(
        self,
        *,
        l: int,
        r: float,
        z: float,
        value: float,
    ):
        self._l = l
        self._r = r
        self._z = z
        self._value = value

    @property
    def l(self):
        return self._l

    @property
    def r(self):
        return self._r

    @property
    def z(self):
        return self._z

    @property
    def value(self):
        return self._value


class AverageMultipoles(Representation):
    def __init__(
        self,
        *,
        l: int,
        r: float,
        z_min: float,
        z_max: float,
        value: float,
    ):
        self._l = l
        self._r = r
        self._z_min = z_min
        self._z_max = z_max
        self._value = value

    @property
    def l(self):
        return self._l

    @property
    def r(self):
        return self._r

    @property
    def z_min(self):
        return self._z_min

    @property
    def z_max(self):
        return self._z_max

    @property
    def value(self):
        return self._value
