import math
import os
import sys


def _grid_tolerance(t, epoch, elapsed, dt):
    # 留浮点容差，上限 dt / 32
    representation = 4 * max(math.ulp(t), math.ulp(epoch), math.ulp(elapsed))
    accumulated = 128 * sys.float_info.epsilon * abs(elapsed)
    return min(dt / 32, max(representation, accumulated))


def validate_integration_request(t, *, epoch, dt, completed_steps=0,
                                 checkpoint=None, checkpoint_interval=None):
    if dt is None:
        raise ValueError("sim.dt is not set")
    dt = float(dt)
    if not math.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be finite and positive")
    t, epoch = float(t), float(epoch)
    if not math.isfinite(t) or not math.isfinite(epoch):
        raise ValueError("target and epoch must be finite")
    # 从初始 epoch 数步，不重新起网格
    steps = (t - epoch) / dt
    if not math.isfinite(steps):
        raise ValueError(f"(target - epoch)/dt must be finite (target={t:.17g}, epoch={epoch:.17g}, dt={dt:.17g})")
    nearest = round(steps)
    elapsed = nearest * dt
    aligned = epoch + elapsed
    # 注意相邻时刻是否还能区分
    if not (epoch + (nearest - 1) * dt < aligned
            < epoch + (nearest + 1) * dt):
        raise ValueError(
            f"dt={dt:.17g} cannot distinguish adjacent times at epoch={epoch:.17g}, target={t:.17g}; "
            "use a relative initial epoch and suitable time units"
        )
    if abs(t - aligned) > _grid_tolerance(t, epoch, elapsed, dt):
        raise ValueError(f"target={t:.17g} is not aligned to dt={dt:.17g}; nearest allowed target={aligned:.17g}")
    if nearest < completed_steps:
        raise ValueError(f"target={t:.17g} must be at or after the current time")
    # Component integrators and split-stage arithmetic still use absolute times.
    # Keep the outer step above the supported resolution at the largest epoch.
    time_resolution = 1e-14 * max(1.0, abs(epoch), abs(t))
    if nearest > completed_steps and dt <= time_resolution:
        raise ValueError(
            f"dt={dt:.17g} must exceed native time resolution={time_resolution:.17g}; "
            "use a relative initial epoch and suitable time units"
        )
    if checkpoint_interval is not None and checkpoint is None:
        raise ValueError("checkpoint_interval requires a checkpoint filename")
    # interval 换成步数，只在完整外层步后保存
    cadence = 1024
    if checkpoint is not None:
        # Check the path protocol now, without creating or replacing a file.
        path = os.fspath(checkpoint)
        if not isinstance(path, str) or not path or "\0" in path:
            raise ValueError("checkpoint must be a nonempty text path without NUL characters")
        if checkpoint_interval is not None:
            interval = float(checkpoint_interval)
            if not math.isfinite(interval) or interval <= 0:
                raise ValueError("checkpoint_interval must be finite and positive")
            ratio = interval / dt
            if not math.isfinite(ratio):
                raise ValueError(f"checkpoint_interval/dt must be finite (checkpoint_interval={interval:.17g}, dt={dt:.17g})")
            cadence = round(ratio)
            aligned_interval = cadence * dt
            if (cadence < 1
                    or not (cadence - 1) * dt < aligned_interval < (cadence + 1) * dt
                    or abs(interval - aligned_interval)
                    > _grid_tolerance(interval, 0, aligned_interval, dt)):
                raise ValueError(f"checkpoint_interval={interval:.17g} must be a resolvable positive integer multiple of dt={dt:.17g}")
    return nearest, cadence
