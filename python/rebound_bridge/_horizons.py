"""Model-local Horizons context around the official REBOUND query path."""
# 复用 REBOUND 的星历，统一日期和参考平面
import datetime
import re
import threading

import rebound


_QUERY_LOCK = threading.Lock()


def _date_string(value):
    if isinstance(value, datetime.datetime):
        if value.utcoffset() is not None:
            raise ValueError("Horizons datetime date must be timezone-naive")
        return value.strftime("%Y-%m-%d %H:%M:%S")
    if not isinstance(value, str):
        raise TypeError("Horizons date must be a string or timezone-naive datetime")
    if value.startswith("JD"):
        if not re.fullmatch(r"JD[0-9]+(?:\.[0-9]*)?", value):
            raise ValueError("Horizons Julian date must use JD followed by an unsigned decimal number")
        return value
    for pattern in ("%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.datetime.strptime(value, pattern).strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            pass
    raise ValueError("Horizons date must use YYYY-MM-DD[ HH:MM[:SS]] or JD followed by an unsigned decimal number")


def validate_context(context):
    """Validate optional checkpoint provenance without making a remote query."""
    if not isinstance(context, dict) or set(context) != {"date", "plane"}:
        raise ValueError("Horizons context must be a dict with exactly 'date' and 'plane'")
    if context["plane"] not in ("ecliptic", "frame"):
        raise ValueError("Horizons plane must be 'ecliptic' or 'frame'")
    try:
        date = _date_string(context["date"])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid Horizons context date: {exc}") from exc
    return {"date": date, "plane": context["plane"]}


def query_particle(source, context, target, *, date=None, plane=None, m=None):

    if not all(source.units.values()):
        raise ValueError("Horizons queries require length, time, and mass units")
    if not isinstance(target, str) or not target.strip():
        raise ValueError("Horizons target must be a nonempty string")
    if context is None:
        if date is None:
            raise ValueError("first Horizons body requires an explicit date")
        normalized_date = _date_string(date)
        resolved_plane = "ecliptic" if plane is None else plane
    else:
        context = validate_context(context)
        previous_date = context["date"]
        normalized_date = previous_date if date is None else _date_string(date)
        resolved_plane = context["plane"] if plane is None else plane
        if normalized_date != previous_date:
            raise ValueError(f"Horizons date {normalized_date!r} must match model date {previous_date!r}")
        if resolved_plane != context["plane"]:
            raise ValueError(f"Horizons plane {resolved_plane!r} must match model plane {context['plane']!r}")
    if resolved_plane not in ("ecliptic", "frame"):
        raise ValueError("Horizons plane must be 'ecliptic' or 'frame'")
    temporary = rebound.Simulation()
    temporary.units = source.units
    options = {"date": normalized_date, "plane": resolved_plane}
    if m is not None:
        options["m"] = m
    with _QUERY_LOCK:
        previous_default = rebound.horizons.INITDATE
        try:
            temporary.add(target, **options)
        finally:
            rebound.horizons.INITDATE = previous_default
    if temporary.N != 1:
        raise ValueError(f"Horizons target must produce exactly one particle (received {temporary.N})")
    particle = temporary.particles[0].copy()
    # 注意清掉临时模拟的指针
    particle._sim = None
    return particle, {"date": normalized_date, "plane": resolved_plane}
