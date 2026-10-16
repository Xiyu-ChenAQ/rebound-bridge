"""Build a physical tree, then compile it into a fixed Bridge hierarchy."""

# 先建物理树，build 时再拆分

import ctypes
import math
from types import MappingProxyType

import rebound

from .bridge import Bridge
from .checkpoint import CheckpointError, _TREE_SCHEMA_VERSION
from .hierarchy import _COORDINATES, _ORBIT_ARGUMENTS as _ORBITAL, checked_integrator, positive
from ._integration import validate_integration_request


class Body:
    """A physical body. ``body.add(...)`` places a new member around this body."""

    # Body 不另存坐标
    # 限定属性，避免拼写错误
    __slots__ = ("_owner", "_name", "_index", "_hash", "_parent", "_children", "_dt", "_integrator")

    def __init__(self, owner, name, index, particle_hash, parent=None):
        self._owner = owner
        self._name = name
        self._index = index
        self._hash = particle_hash
        self._parent = parent
        self._children = []
        self._dt = None
        self._integrator = None

    @property
    def name(self):
        return self._name

    @property
    def parent(self):
        return self._parent

    @property
    def children(self):
        return tuple(self._children)

    @property
    def m(self):
        """Physical mass; editable before build without recomputing orbits."""
        self._owner._require_open()
        # 注意取真实天体质量，不能取代理总质量
        source = self._owner._component(self)
        # c_uint 按 hash 查，int 按下标查
        return source.particles[ctypes.c_uint(self._hash)].m

    @m.setter
    def m(self, value):
        self._owner._require_editable()
        value = float(value)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"Body {self.name!r} mass must be finite and non-negative")
        if self._children and value == 0:
            raise ValueError(f"Parent {self.name!r} must have positive mass")
        self._owner._source.particles[ctypes.c_uint(self._hash)].m = value

    @property
    def dt(self):
        """Internal step for this body's resolved system; None inherits dt_sub."""
        return self._dt

    @dt.setter
    def dt(self, value):
        self._owner._require_editable()
        self._dt = None if value is None else positive(value, f"{self.name}.dt")

    @property
    def integrator(self):
        return self._integrator or self._owner.integrator

    @integrator.setter
    def integrator(self, value):
        self._owner._require_editable()
        self._integrator = self._owner._check_integrator(value)

    @property
    def effective_integrator(self):
        """Integrator of the component that actually evolves this body."""
        self._owner._require_open()
        if self._owner._bridge is not None:
            return self._owner._component(self).integrator
        host = self._owner._enclosing_host(self)
        return host.integrator if host is not None else self._owner.integrator

    def add(self, name=None, **particle):
        """Add a child with REBOUND elements or parent-relative Cartesian state."""
        return self._owner._add(self, name, particle)

    def add_horizons(self, target, *, date=None, name=None, plane=None, m=None):
        """Query a child through REBOUND's Horizons interface.

        The parent sets hierarchy membership; Horizons already supplies global
        coordinates. Date and reference plane inherit the first query's values.
        """
        return self._owner._add_horizons(self, target, date=date, name=name, plane=plane, m=m)

    def __repr__(self):
        return f"Body(name={self.name!r}, children={len(self.children)})"


class Simulation:
    """A physical hierarchy with REBOUND-style particle arguments.

    Use ``sim.add`` for a root body and ``parent.add`` for its children.
    Set ``dt`` for outer coupling and ``dt_sub`` for resolved systems; a host's
    ``dt`` overrides the latter. All descendants of a root's direct child share
    one internal simulation. The physical tree does not introduce additional
    recursive Bridge time levels.
    """

    __slots__ = (
        "_bridge", "_closed", "_source", "_bodies", "_roots", "_hashes", "_dt",
        "_dt_sub", "_integrator", "_method", "_coupling", "_horizons_context",
    )

    def __init__(self, *, units=None, dt=None, dt_sub=None, integrator="whfast", method="lf4_2", coupling="monopole"):
        self._bridge = None
        self._closed = False
        # 初始状态统一存全局坐标
        # build 后读 _bridge，_source 不再更新
        self._source = rebound.Simulation()
        self._bodies = {}
        self._roots = []
        self._hashes = set()
        self._horizons_context = None
        self._dt = self._dt_sub = None
        self._integrator = "whfast"
        self.dt = dt
        self.dt_sub = dt_sub
        self.integrator = integrator
        self.method = method
        self.coupling = coupling
        if units is not None:
            self.units = units

    def _require_open(self):
        if self._closed:
            raise RuntimeError("Simulation is closed")

    def _require_editable(self):
        self._require_open()
        if self._bridge is not None:
            raise RuntimeError("Hierarchy and integration settings are fixed after build(), integrate(), or save_to_file()")

    @staticmethod
    def _check_integrator(value):
        return checked_integrator(value)

    @property
    def G(self):
        """Gravitational constant; changing it preserves existing Cartesian states."""
        self._require_open()
        return self._source.G if self._bridge is None else self._bridge.main.G

    @G.setter
    def G(self, value):
        self._require_editable()
        self._source.G = positive(value, "G")

    @property
    def units(self):
        return self._source.units

    @units.setter
    def units(self, value):
        self._require_editable()
        if self._bodies:
            raise ValueError("units must be set before adding bodies")
        self._source.units = value

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, value):
        self._require_editable()
        self._dt = None if value is None else positive(value, "dt")

    @property
    def dt_sub(self):
        return self._dt_sub

    @dt_sub.setter
    def dt_sub(self, value):
        self._require_editable()
        self._dt_sub = None if value is None else positive(value, "dt_sub")

    @property
    def integrator(self):
        return self._integrator

    @integrator.setter
    def integrator(self, value):
        self._require_editable()
        self._integrator = self._check_integrator(value)

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, value):
        self._require_editable()
        if value not in ("kdk", "lf4_2"):
            raise ValueError("method must be 'kdk' or 'lf4_2'")
        self._method = value

    @property
    def coupling(self):
        return self._coupling

    @coupling.setter
    def coupling(self, value):
        self._require_editable()
        if value not in ("monopole", "resolved"):
            raise ValueError("coupling must be 'monopole' or 'resolved'")
        self._coupling = value

    @property
    def t(self):
        return self._bridge.t if self._bridge is not None else self._source.t

    @t.setter
    def t(self, value):
        self._require_editable()
        value = float(value)
        if not math.isfinite(value):
            raise ValueError("Initial time t must be finite")
        self._source.t = value

    @property
    def bodies(self):
        return MappingProxyType(self._bodies)

    @property
    def N(self):
        """Number of physical bodies, excluding centre-of-mass proxies."""
        return len(self._bodies)

    @property
    def horizons(self):
        """Read-only initial Horizons date/plane, or None before a query."""
        return None if self._horizons_context is None else MappingProxyType(self._horizons_context)

    @property
    def children(self):
        return tuple(self._roots)

    def __getitem__(self, name):
        return self._bodies[name]

    def add(self, name=None, **particle):
        """Add a root body. Use its returned handle to add orbiting children."""
        return self._add(None, name, particle)

    def add_horizons(self, target, *, date=None, name=None, plane=None, m=None):
        """Add a root using REBOUND's Horizons query and unit conversion.

        Set physical units first and give an explicit date on the first query.
        Later queries inherit that date and plane. ``m``, when supplied, uses
        this simulation's mass unit. The initial simulation time stays unchanged.
        """
        return self._add_horizons(None, target, date=date, name=name, plane=plane, m=m)

    def _identity(self, name, supplied_hash=None):
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Body name must be a nonempty string")
        if name in self._bodies:
            raise ValueError(f"Body name {name!r} is already in the hierarchy")
        particle_hash = rebound.hash(name).value
        if supplied_hash is not None and rebound.hash(supplied_hash).value != particle_hash:
            raise ValueError("name and hash must identify the same body; omit hash to use name")
        if not particle_hash or particle_hash in self._hashes:
            raise ValueError(f"Body {name!r} has a zero or duplicate REBOUND hash; choose a different name")
        return particle_hash

    def _parent_particle(self, parent):
        if parent is None:
            return None
        primary = self._source.particles[parent._index]
        if primary.m <= 0:
            raise ValueError(f"Parent {parent.name!r} must have positive mass")
        return primary

    def _add_horizons(self, parent, target, *, date, name, plane, m):
        from ._horizons import query_particle

        self._require_editable()
        name = target if name is None else name
        particle_hash = self._identity(name)
        self._parent_particle(parent)
        candidate, context = query_particle(self._source, self._horizons_context, target,
                                            date=date, plane=plane, m=m)
        candidate.hash = name
        # 已是全局状态，不再加宿主偏移
        body = self._append(parent, name, particle_hash, candidate)
        self._horizons_context = context
        return body

    def _add(self, parent, name, kwargs):
        self._require_editable()
        if "primary" in kwargs:
            raise ValueError("primary must not be supplied; use parent.add(...)")
        if kwargs.get("jacobi_masses", False):
            raise ValueError("jacobi_masses must be false")
        if {"particle", "simulation", "variation", "variation2"}.intersection(kwargs):
            raise ValueError("particle, simulation, variation, and variation2 are not supported; use Bridge.from_simulation() for existing particle data")
        if name is None:
            name = kwargs.get("hash", f"body_{len(self._bodies)}")
        particle_hash = self._identity(name, kwargs.get("hash", name))
        kwargs = {**kwargs, "hash": name}
        orbital = any(kwargs.get(key) is not None for key in _ORBITAL)
        if parent is None and orbital:
            raise ValueError("Orbital elements require a parent; use parent.add(...)")
        if parent is not None:
            primary = self._parent_particle(parent)
            if orbital:
                # 指定宿主，根数转换用 REBOUND
                kwargs["primary"] = primary
            else:
                # 相对宿主的坐标转回全局
                kwargs.update({key: float(0.0 if kwargs.get(key) is None else kwargs[key])
                               + getattr(primary, key) for key in _COORDINATES})
        candidate = rebound.Particle(simulation=self._source, **kwargs)
        return self._append(parent, name, particle_hash, candidate)

    def _append(self, parent, name, particle_hash, candidate):
        if not math.isfinite(candidate.m) or candidate.m < 0:
            raise ValueError(f"Body {name!r} mass must be finite and non-negative")
        if not all(math.isfinite(getattr(candidate, key)) for key in _COORDINATES):
            raise ValueError(f"Body {name!r} position and velocity must be finite")
        body = Body(self, name, self._source.N, particle_hash, parent)
        self._source.add(candidate)
        self._bodies[name] = body
        self._hashes.add(particle_hash)
        (self._roots if parent is None else parent._children).append(body)
        return body

    @staticmethod
    def _enclosing_host(body):
        # 向上找子系统宿主，与 _groups 对应
        while body.parent is not None and body.parent.parent is not None:
            body = body.parent
        return body if body.parent is not None and body.children else None

    def _component(self, body):
        self._require_open()
        if self._bridge is None:
            return self._source
        self._bridge._require_open()
        host = self._enclosing_host(body)
        return self._bridge.main if host is None else self._bridge.subsystems[host.name].simulation

    @staticmethod
    def _descendants(body):
        yield body
        for child in body.children:
            yield from Simulation._descendants(child)

    def _groups(self):
        # 宿主及全部后代分到一起，不再嵌套 Bridge
        return {child.name: tuple(self._descendants(child)) for root in self._roots for child in root.children if child.children}

    def _steps(self, groups):
        if self.dt is None:
            raise ValueError("sim.dt must be set before build() or integrate()")
        steps = {}
        for name, members in groups.items():
            step = members[0].dt if members[0].dt is not None else self.dt_sub
            if step is None:
                raise ValueError(f"Subsystem {name!r} requires sim.dt_sub or sim[{name!r}].dt")
            ratio = self.dt / step
            if not math.isfinite(ratio) or ratio < 1 or abs(ratio - round(ratio)) > 1e-12:
                raise ValueError(f"sim.dt must be an integer multiple of the internal step for {name!r}; got {self.dt!r}/{step!r} = {ratio!r}")
            steps[name] = step
        for body in self._bodies.values():
            if body.name not in groups and (body.dt is not None or body._integrator is not None):
                raise ValueError(f"{body.name!r} is not a subsystem host; set dt and integrator on its subsystem host or Simulation")
        return steps

    def _tree_metadata(self):
        return {
            "schema": _TREE_SCHEMA_VERSION, "dt_sub": self.dt_sub, "integrator": self.integrator,
            "horizons": None if self._horizons_context is None else dict(self._horizons_context),
            "bodies": [{"name": b.name, "hash": b._hash, "parent": b.parent.name if b.parent else None,
                        "dt": b.dt, "integrator": b._integrator} for b in self._bodies.values()],
        }

    def build(self):
        """Compile the physical tree to COM proxies and fixed internal systems."""
        self._require_open()
        if self._bridge is None:
            groups = self._groups()
            steps = self._steps(groups)
            source = self._source.copy()
            source.integrator = self.integrator
            source.dt = self.dt
            # 整理成员索引，构建复用 Bridge
            bridge = Bridge.from_simulation(
                source, subsystems={name: [b._index for b in members] for name, members in groups.items()},
                dt=self.dt, dt_sub=steps, method=self.method, coupling=self.coupling,
                integrator=self.integrator,
            )
            try:
                for name, members in groups.items():
                    bridge.subsystems[name].integrator = members[0].integrator
                bridge.build()
                bridge._tree_metadata = self._tree_metadata()
            except BaseException:
                bridge.close()
                raise
            # 成功后再接管
            self._bridge = bridge
        else:
            self._bridge.build()
        return self

    def integrate(self, t, *, checkpoint=None, checkpoint_interval=None):
        """Advance to absolute simulation time ``t`` on the outer-step grid.

        The first call builds and freezes the hierarchy. ``t`` must lie on the
        grid defined by the initial time and ``dt``. Optional checkpoints save
        completed outer steps; ``checkpoint_interval`` uses simulation time
        units. After a failed integration, restore the last saved checkpoint.
        """
        self._require_open()
        if self._bridge is None:
            if self.dt is None:
                raise ValueError("sim.dt must be set before build() or integrate()")
            # 先检查时间，别提前冻结层级
            validate_integration_request(
                t, epoch=self._source.t, dt=self.dt,
                checkpoint=checkpoint, checkpoint_interval=checkpoint_interval,
            )
        self.build()
        self._bridge.integrate(t, checkpoint=checkpoint, checkpoint_interval=checkpoint_interval)

    def to_rebound(self):
        """Return a detached physical snapshot without freezing an unbuilt tree."""
        self._require_open()
        return self._bridge.to_rebound() if self._bridge is not None else self._source.copy()

    def save_to_file(self, filename):
        self.build()
        self._bridge.save_to_file(filename)

    @classmethod
    def from_file(cls, filename):
        bridge = Bridge.from_file(filename)
        try:
            tree = bridge._tree_metadata
            if not isinstance(tree, dict) or not isinstance(tree.get("bodies"), list):
                raise ValueError("Checkpoint lacks a physical tree; use Bridge.from_file()")
            sim = cls(dt=bridge.dt, dt_sub=tree["dt_sub"], integrator=tree["integrator"], method=bridge.method, coupling=bridge.coupling)
            from ._horizons import validate_context
            context = tree.get("horizons")
            if context is not None:
                sim._horizons_context = validate_context(context)
            sim._source = bridge.to_rebound()
            # 只恢复树结构，不重新初始化积分器
            for record in tree["bodies"]:
                name, key, parent_name = record["name"], record["hash"], record["parent"]
                if not isinstance(name, str) or not name.strip() or name in sim._bodies:
                    raise ValueError(f"Invalid or duplicate body name {name!r} in physical tree")
                if type(key) is not int or key != rebound.hash(name).value or not key or key in sim._hashes:
                    raise ValueError(f"Invalid or duplicate particle hash {key!r} for {name!r} in physical tree")
                parent = None if parent_name is None else sim._bodies[parent_name]
                index = sim._source.particles[ctypes.c_uint(key)].index
                body = Body(sim, name, index, key, parent)
                body.dt = record["dt"]
                if record["integrator"] is not None:
                    body.integrator = record["integrator"]
                sim._bodies[name] = body
                sim._hashes.add(key)
                (sim._roots if parent is None else parent._children).append(body)
            physical_hashes = [p.hash.value for p in sim._source.particles]
            if len(physical_hashes) != len(sim._bodies) or set(physical_hashes) != sim._hashes:
                raise ValueError("Physical tree does not match the stored particles")
            groups = sim._groups()
            steps = sim._steps(groups)
            if set(groups) != set(bridge.subsystems) or sim.integrator != bridge.main.integrator:
                raise ValueError("Physical tree does not match the saved integration components")
            for name, members in groups.items():
                sub = bridge.subsystems[name]
                hashes = [p.hash.value for p in sub.particles]
                if (set(hashes) != {b._hash for b in members} or hashes[0] != members[0]._hash
                        or sub.dt != steps[name] or sub.integrator != members[0].integrator
                        or bridge.main.particles[sub.host_index].hash.value != members[0]._hash):
                    raise ValueError(f"Physical tree and saved subsystem {name!r} disagree")
            sim._bridge = bridge
            return sim
        except Exception as exc:
            bridge.close()
            raise CheckpointError(f"Invalid physical hierarchy checkpoint: {exc}") from exc

    def hierarchy(self):
        """Return a readable tree of the physical body names."""
        lines = []

        def visit(body, prefix):
            for i, child in enumerate(body.children):
                last = i == len(body.children) - 1
                lines.append(prefix + ("└── " if last else "├── ") + child.name)
                visit(child, prefix + ("    " if last else "│   "))

        for root in self._roots:
            lines.append(root.name)
            visit(root, "")
        return "\n".join(lines)

    def status(self):
        """Print time, physical hierarchy, and component integration settings."""
        self._require_open()
        if self._bridge is not None:
            self._bridge._require_open()
        print(f"rebound_bridge Simulation: {self.N} physical bodies, t={self.t:g}")
        print(f"Method: {self.method}; coupling: {self.coupling}; outer dt={self.dt}")
        print(f"Main integrator: {self.integrator}")
        for name, members in self._groups().items():
            host = members[0]
            step = self.dt_sub if host.dt is None else host.dt
            print(f"Subsystem {name}: {len(members)} bodies, {host.effective_integrator}, dt={step}")
        print(self.hierarchy())

    def close(self):
        if self._bridge is not None:
            self._bridge.close()
        self._closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()

    def __repr__(self):
        return f"Simulation(t={self.t:g}, bodies={len(self._bodies)}, dt={self.dt!r})"
