"""Minimal Surface description for the certified layer (extended in a later task)."""


class Surface:
    def __init__(self, vertices_hyp, is_bolza=False, label="custom"):
        self.vertices_hyp = vertices_hyp
        self.is_bolza = is_bolza
        self.label = label

    @classmethod
    def bolza(cls) -> "Surface":
        from ._arb import prec
        from .mesh import bolza_vertices_hyp
        # Vertex balls must be built at the certification precision (96 bits); default-precision balls
        # (53 bits) would cap the whole certified computation at ~1e-8 absolute radius.
        with prec(96):
            return cls(bolza_vertices_hyp(), is_bolza=True, label="bolza")
