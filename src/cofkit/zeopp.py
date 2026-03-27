from __future__ import annotations

import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

COFKIT_ZEOPP_ENV_VAR = "COFKIT_ZEOPP_PATH"


class ZeoppError(RuntimeError):
    """Base error for Zeo++ wrapper failures."""


class ZeoppConfigurationError(ZeoppError):
    """Raised when the Zeo++ binary cannot be resolved."""


class ZeoppExecutionError(ZeoppError):
    """Raised when a Zeo++ command exits unsuccessfully."""


class ZeoppParseError(ZeoppError):
    """Raised when a Zeo++ output file cannot be parsed."""


@dataclass(frozen=True)
class ZeoppBasicPoreProperties:
    largest_included_sphere: float
    largest_free_sphere: float
    largest_included_sphere_along_free_path: float
    axis_aligned_free_sphere: Mapping[str, float]
    axis_aligned_included_sphere_along_free_path: Mapping[str, float]

    def to_dict(self) -> dict[str, object]:
        return {
            "largest_included_sphere": self.largest_included_sphere,
            "largest_free_sphere": self.largest_free_sphere,
            "largest_included_sphere_along_free_path": self.largest_included_sphere_along_free_path,
            "axis_aligned_free_sphere": dict(self.axis_aligned_free_sphere),
            "axis_aligned_included_sphere_along_free_path": dict(
                self.axis_aligned_included_sphere_along_free_path
            ),
        }


@dataclass(frozen=True)
class ZeoppChannelEntry:
    index: int
    largest_included_sphere: float
    largest_free_sphere: float
    largest_included_sphere_along_free_path: float

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "largest_included_sphere": self.largest_included_sphere,
            "largest_free_sphere": self.largest_free_sphere,
            "largest_included_sphere_along_free_path": self.largest_included_sphere_along_free_path,
        }


@dataclass(frozen=True)
class ZeoppChannelSummary:
    n_channels: int
    n_pockets: int | None
    channel_dimensionality: int | None
    channel_dimensionalities: tuple[int, ...]
    probe_radius: float
    probe_diameter: float | None
    largest_included_sphere: float | None
    largest_free_sphere: float | None
    largest_included_sphere_along_free_path: float | None
    channels: tuple[ZeoppChannelEntry, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "n_channels": self.n_channels,
            "n_pockets": self.n_pockets,
            "channel_dimensionality": self.channel_dimensionality,
            "channel_dimensionalities": list(self.channel_dimensionalities),
            "probe_radius": self.probe_radius,
            "probe_diameter": self.probe_diameter,
            "largest_included_sphere": self.largest_included_sphere,
            "largest_free_sphere": self.largest_free_sphere,
            "largest_included_sphere_along_free_path": self.largest_included_sphere_along_free_path,
            "channels": [channel.to_dict() for channel in self.channels],
        }


@dataclass(frozen=True)
class ZeoppSurfaceAreaProperties:
    unitcell_volume: float
    density: float
    accessible_surface_area_a2: float
    accessible_surface_area_m2_cm3: float
    accessible_surface_area_m2_g: float
    non_accessible_surface_area_a2: float
    non_accessible_surface_area_m2_cm3: float
    non_accessible_surface_area_m2_g: float
    n_channels: int
    channel_surface_area_a2: tuple[float, ...]
    n_pockets: int
    pocket_surface_area_a2: tuple[float, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "unitcell_volume": self.unitcell_volume,
            "density": self.density,
            "accessible_surface_area_a2": self.accessible_surface_area_a2,
            "accessible_surface_area_m2_cm3": self.accessible_surface_area_m2_cm3,
            "accessible_surface_area_m2_g": self.accessible_surface_area_m2_g,
            "non_accessible_surface_area_a2": self.non_accessible_surface_area_a2,
            "non_accessible_surface_area_m2_cm3": self.non_accessible_surface_area_m2_cm3,
            "non_accessible_surface_area_m2_g": self.non_accessible_surface_area_m2_g,
            "n_channels": self.n_channels,
            "channel_surface_area_a2": list(self.channel_surface_area_a2),
            "n_pockets": self.n_pockets,
            "pocket_surface_area_a2": list(self.pocket_surface_area_a2),
        }


@dataclass(frozen=True)
class ZeoppVolumeProperties:
    unitcell_volume: float
    density: float
    accessible_volume_a3: float
    accessible_volume_fraction: float
    accessible_volume_cm3_g: float
    non_accessible_volume_a3: float
    non_accessible_volume_fraction: float
    non_accessible_volume_cm3_g: float
    n_channels: int
    channel_volume_a3: tuple[float, ...]
    n_pockets: int
    pocket_volume_a3: tuple[float, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "unitcell_volume": self.unitcell_volume,
            "density": self.density,
            "accessible_volume_a3": self.accessible_volume_a3,
            "accessible_volume_fraction": self.accessible_volume_fraction,
            "accessible_volume_cm3_g": self.accessible_volume_cm3_g,
            "non_accessible_volume_a3": self.non_accessible_volume_a3,
            "non_accessible_volume_fraction": self.non_accessible_volume_fraction,
            "non_accessible_volume_cm3_g": self.non_accessible_volume_cm3_g,
            "n_channels": self.n_channels,
            "channel_volume_a3": list(self.channel_volume_a3),
            "n_pockets": self.n_pockets,
            "pocket_volume_a3": list(self.pocket_volume_a3),
        }


@dataclass(frozen=True)
class ZeoppAccessibilitySummary:
    n_voronoi_nodes: int
    n_accessible_nodes: int
    n_inaccessible_nodes: int
    accessible_fraction: float | None

    def to_dict(self) -> dict[str, object]:
        return {
            "n_voronoi_nodes": self.n_voronoi_nodes,
            "n_accessible_nodes": self.n_accessible_nodes,
            "n_inaccessible_nodes": self.n_inaccessible_nodes,
            "accessible_fraction": self.accessible_fraction,
        }


@dataclass(frozen=True)
class ZeoppProbeScanSettings:
    channel_radius: float
    probe_radius: float
    surface_samples_per_atom: int
    volume_samples_total: int

    def to_dict(self) -> dict[str, object]:
        return {
            "channel_radius": self.channel_radius,
            "probe_radius": self.probe_radius,
            "surface_samples_per_atom": self.surface_samples_per_atom,
            "volume_samples_total": self.volume_samples_total,
        }


@dataclass(frozen=True)
class ZeoppProbeScanResult:
    settings: ZeoppProbeScanSettings
    output_dir: str
    output_paths: Mapping[str, str]
    status: str
    error: str | None = None
    channel_summary: ZeoppChannelSummary | None = None
    surface_area: ZeoppSurfaceAreaProperties | None = None
    volume: ZeoppVolumeProperties | None = None
    accessibility: ZeoppAccessibilitySummary | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "settings": self.settings.to_dict(),
            "output_dir": self.output_dir,
            "output_paths": dict(self.output_paths),
            "status": self.status,
            "error": self.error,
            "channel_summary": None if self.channel_summary is None else self.channel_summary.to_dict(),
            "surface_area": None if self.surface_area is None else self.surface_area.to_dict(),
            "volume": None if self.volume is None else self.volume.to_dict(),
            "accessibility": None if self.accessibility is None else self.accessibility.to_dict(),
        }


@dataclass(frozen=True)
class ZeoppBaselineResult:
    basic_pore_properties: ZeoppBasicPoreProperties
    point_probe_channels: ZeoppChannelSummary
    point_probe_surface_area: ZeoppSurfaceAreaProperties
    point_probe_volume: ZeoppVolumeProperties
    output_paths: Mapping[str, str]

    def to_dict(self) -> dict[str, object]:
        return {
            "basic_pore_properties": self.basic_pore_properties.to_dict(),
            "point_probe_channels": self.point_probe_channels.to_dict(),
            "point_probe_surface_area": self.point_probe_surface_area.to_dict(),
            "point_probe_volume": self.point_probe_volume.to_dict(),
            "output_paths": dict(self.output_paths),
        }


@dataclass(frozen=True)
class ZeoppAnalysisResult:
    input_cif: str
    zeopp_binary: str
    output_dir: str
    report_path: str
    baseline: ZeoppBaselineResult
    probe_scans: tuple[ZeoppProbeScanResult, ...] = ()

    @property
    def properties(self) -> ZeoppBasicPoreProperties:
        return self.baseline.basic_pore_properties

    def to_dict(self) -> dict[str, object]:
        return {
            "input_cif": self.input_cif,
            "zeopp_binary": self.zeopp_binary,
            "output_dir": self.output_dir,
            "report_path": self.report_path,
            "baseline": self.baseline.to_dict(),
            "probe_scans": [scan.to_dict() for scan in self.probe_scans],
        }


def resolve_zeopp_binary(zeopp_path: str | Path | None = None) -> Path:
    raw_value = str(zeopp_path) if zeopp_path is not None else os.environ.get(COFKIT_ZEOPP_ENV_VAR)
    if not raw_value:
        raise ZeoppConfigurationError(
            f"Zeo++ binary is not configured. Set {COFKIT_ZEOPP_ENV_VAR} to the Zeo++ network binary path."
        )

    path = Path(raw_value).expanduser()
    if not path.is_file():
        raise ZeoppConfigurationError(f"Zeo++ binary does not exist: {path}")
    if not os.access(path, os.X_OK):
        raise ZeoppConfigurationError(f"Zeo++ binary is not executable: {path}")
    return path.resolve()


def analyze_zeopp_pore_properties(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    probe_radii: Sequence[float] = (),
    channel_radius: float | None = None,
    surface_samples_per_atom: int = 250,
    volume_samples_total: int = 5000,
    zeopp_path: str | Path | None = None,
    timeout_seconds: float = 300.0,
    continue_on_probe_error: bool = True,
) -> ZeoppAnalysisResult:
    input_path = Path(cif_path).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"CIF file does not exist: {input_path}")
    if timeout_seconds <= 0.0:
        raise ValueError("timeout_seconds must be positive.")
    if surface_samples_per_atom <= 0:
        raise ValueError("surface_samples_per_atom must be positive.")
    if volume_samples_total <= 0:
        raise ValueError("volume_samples_total must be positive.")

    normalized_probe_radii = tuple(float(value) for value in probe_radii)
    for value in normalized_probe_radii:
        if value < 0.0:
            raise ValueError("probe_radii values must be non-negative.")
    if channel_radius is not None and channel_radius < 0.0:
        raise ValueError("channel_radius must be non-negative when provided.")

    binary = resolve_zeopp_binary(zeopp_path)
    run_dir = _resolve_output_dir(input_path, output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    baseline = _run_baseline_analysis(
        binary,
        input_path,
        run_dir / "baseline",
        surface_samples_per_atom=surface_samples_per_atom,
        volume_samples_total=volume_samples_total,
        timeout_seconds=timeout_seconds,
    )

    probe_scans = tuple(
        _run_probe_scan(
            binary,
            input_path,
            run_dir / "probe_scans",
            scan_index=index,
            settings=ZeoppProbeScanSettings(
                channel_radius=value if channel_radius is None else float(channel_radius),
                probe_radius=value,
                surface_samples_per_atom=surface_samples_per_atom,
                volume_samples_total=volume_samples_total,
            ),
            timeout_seconds=timeout_seconds,
            continue_on_error=continue_on_probe_error,
        )
        for index, value in enumerate(normalized_probe_radii, start=1)
    )

    result = ZeoppAnalysisResult(
        input_cif=str(input_path),
        zeopp_binary=str(binary),
        output_dir=str(run_dir),
        report_path=str(run_dir / "zeopp_report.json"),
        baseline=baseline,
        probe_scans=probe_scans,
    )
    Path(result.report_path).write_text(json.dumps(result.to_dict(), indent=2), encoding="utf-8")
    return result


def analyze_zeopp_basic_pore_properties(
    cif_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    probe_radius: float = 1.86,
    zeopp_path: str | Path | None = None,
    timeout_seconds: float = 300.0,
) -> ZeoppAnalysisResult:
    return analyze_zeopp_pore_properties(
        cif_path,
        output_dir=output_dir,
        probe_radii=(probe_radius,),
        zeopp_path=zeopp_path,
        timeout_seconds=timeout_seconds,
    )


def _run_baseline_analysis(
    binary: Path,
    input_path: Path,
    output_dir: Path,
    *,
    surface_samples_per_atom: int,
    volume_samples_total: int,
    timeout_seconds: float,
) -> ZeoppBaselineResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    output_paths = _command_output_paths(
        output_dir,
        {
            "res": f"{stem}.res",
            "resex": f"{stem}.resex",
            "chan": f"{stem}.chan",
            "sa": f"{stem}.sa",
            "vol": f"{stem}.vol",
        },
    )

    _run_zeopp_command(
        binary,
        ("-res", output_paths["res"], str(input_path)),
        stdout_log_path=Path(output_paths["res_stdout_log"]),
        stderr_log_path=Path(output_paths["res_stderr_log"]),
        timeout_seconds=timeout_seconds,
    )
    _run_zeopp_command(
        binary,
        ("-resex", output_paths["resex"], str(input_path)),
        stdout_log_path=Path(output_paths["resex_stdout_log"]),
        stderr_log_path=Path(output_paths["resex_stderr_log"]),
        timeout_seconds=timeout_seconds,
    )
    chan_stdout = _run_zeopp_command(
        binary,
        ("-chan", "0", output_paths["chan"], str(input_path)),
        stdout_log_path=Path(output_paths["chan_stdout_log"]),
        stderr_log_path=Path(output_paths["chan_stderr_log"]),
        timeout_seconds=timeout_seconds,
    )
    _run_zeopp_command(
        binary,
        ("-sa", "0", "0", str(surface_samples_per_atom), output_paths["sa"], str(input_path)),
        stdout_log_path=Path(output_paths["sa_stdout_log"]),
        stderr_log_path=Path(output_paths["sa_stderr_log"]),
        timeout_seconds=timeout_seconds,
    )
    _run_zeopp_command(
        binary,
        ("-vol", "0", "0", str(volume_samples_total), output_paths["vol"], str(input_path)),
        stdout_log_path=Path(output_paths["vol_stdout_log"]),
        stderr_log_path=Path(output_paths["vol_stderr_log"]),
        timeout_seconds=timeout_seconds,
    )

    return ZeoppBaselineResult(
        basic_pore_properties=_parse_resex_output(Path(output_paths["resex"])),
        point_probe_channels=_parse_chan_output(Path(output_paths["chan"]), chan_stdout),
        point_probe_surface_area=_parse_surface_area_output(Path(output_paths["sa"])),
        point_probe_volume=_parse_volume_output(Path(output_paths["vol"])),
        output_paths=output_paths,
    )


def _run_probe_scan(
    binary: Path,
    input_path: Path,
    output_root: Path,
    *,
    scan_index: int,
    settings: ZeoppProbeScanSettings,
    timeout_seconds: float,
    continue_on_error: bool,
) -> ZeoppProbeScanResult:
    label = (
        f"probe_scan_{scan_index:02d}"
        f"__chan_{_format_radius_label(settings.channel_radius)}"
        f"__probe_{_format_radius_label(settings.probe_radius)}"
    )
    output_dir = output_root / label
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    output_paths = _command_output_paths(
        output_dir,
        {
            "chan": f"{stem}.chan",
            "sa": f"{stem}.sa",
            "vol": f"{stem}.vol",
            "axs": f"{stem}.axs",
        },
    )

    errors: list[str] = []
    channel_summary: ZeoppChannelSummary | None = None
    surface_area: ZeoppSurfaceAreaProperties | None = None
    volume: ZeoppVolumeProperties | None = None
    accessibility: ZeoppAccessibilitySummary | None = None

    try:
        chan_stdout = _run_zeopp_command(
            binary,
            ("-chan", f"{settings.probe_radius:g}", output_paths["chan"], str(input_path)),
            stdout_log_path=Path(output_paths["chan_stdout_log"]),
            stderr_log_path=Path(output_paths["chan_stderr_log"]),
            timeout_seconds=timeout_seconds,
        )
        channel_summary = _parse_chan_output(Path(output_paths["chan"]), chan_stdout)
    except ZeoppError as exc:
        errors.append(f"chan: {exc}")
        if not continue_on_error:
            raise

    try:
        _run_zeopp_command(
            binary,
            (
                "-sa",
                f"{settings.channel_radius:g}",
                f"{settings.probe_radius:g}",
                str(settings.surface_samples_per_atom),
                output_paths["sa"],
                str(input_path),
            ),
            stdout_log_path=Path(output_paths["sa_stdout_log"]),
            stderr_log_path=Path(output_paths["sa_stderr_log"]),
            timeout_seconds=timeout_seconds,
        )
        surface_area = _parse_surface_area_output(Path(output_paths["sa"]))
    except ZeoppError as exc:
        errors.append(f"sa: {exc}")
        if not continue_on_error:
            raise

    try:
        _run_zeopp_command(
            binary,
            (
                "-vol",
                f"{settings.channel_radius:g}",
                f"{settings.probe_radius:g}",
                str(settings.volume_samples_total),
                output_paths["vol"],
                str(input_path),
            ),
            stdout_log_path=Path(output_paths["vol_stdout_log"]),
            stderr_log_path=Path(output_paths["vol_stderr_log"]),
            timeout_seconds=timeout_seconds,
        )
        volume = _parse_volume_output(Path(output_paths["vol"]))
    except ZeoppError as exc:
        errors.append(f"vol: {exc}")
        if not continue_on_error:
            raise

    try:
        _run_zeopp_command(
            binary,
            ("-axs", f"{settings.probe_radius:g}", output_paths["axs"], str(input_path)),
            stdout_log_path=Path(output_paths["axs_stdout_log"]),
            stderr_log_path=Path(output_paths["axs_stderr_log"]),
            timeout_seconds=timeout_seconds,
        )
        accessibility = _parse_accessibility_output(Path(output_paths["axs"]))
    except ZeoppError as exc:
        errors.append(f"axs: {exc}")
        if not continue_on_error:
            raise

    return ZeoppProbeScanResult(
        settings=settings,
        output_dir=str(output_dir),
        output_paths=output_paths,
        status="ok" if not errors else "error",
        error=None if not errors else "; ".join(errors),
        channel_summary=channel_summary,
        surface_area=surface_area,
        volume=volume,
        accessibility=accessibility,
    )


def _resolve_output_dir(input_path: Path, output_dir: str | Path | None) -> Path:
    if output_dir is None:
        return input_path.parent / f"{input_path.stem}_zeopp"
    return Path(output_dir).expanduser().resolve()


def _command_output_paths(output_dir: Path, files: Mapping[str, str]) -> dict[str, str]:
    output_paths: dict[str, str] = {}
    for key, filename in files.items():
        output_paths[key] = str(output_dir / filename)
        output_paths[f"{key}_stdout_log"] = str(output_dir / f"zeopp_{key}.stdout.log")
        output_paths[f"{key}_stderr_log"] = str(output_dir / f"zeopp_{key}.stderr.log")
    return output_paths


def _run_zeopp_command(
    binary: Path,
    args: tuple[str, ...],
    *,
    stdout_log_path: Path,
    stderr_log_path: Path,
    timeout_seconds: float,
) -> str:
    command = [str(binary), *args]
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout if isinstance(exc.stdout, str) else (exc.stdout or b"").decode("utf-8", errors="replace")
        stderr = exc.stderr if isinstance(exc.stderr, str) else (exc.stderr or b"").decode("utf-8", errors="replace")
        stdout_log_path.write_text(stdout, encoding="utf-8")
        stderr_log_path.write_text(stderr, encoding="utf-8")
        raise ZeoppExecutionError(
            f"Zeo++ command timed out after {timeout_seconds:g} seconds: {' '.join(command)}"
        ) from exc

    stdout_log_path.write_text(completed.stdout or "", encoding="utf-8")
    stderr_log_path.write_text(completed.stderr or "", encoding="utf-8")
    if completed.returncode != 0:
        raise ZeoppExecutionError(
            "Zeo++ command failed with "
            f"exit code {completed.returncode}: {' '.join(command)}"
        )
    return completed.stdout or ""


def _parse_resex_output(resex_output_path: Path) -> ZeoppBasicPoreProperties:
    if not resex_output_path.is_file():
        raise ZeoppParseError(f"Zeo++ -resex output file was not created: {resex_output_path}")

    for line in reversed(resex_output_path.read_text(encoding="utf-8").splitlines()):
        floats = _parse_trailing_floats(line)
        if len(floats) >= 9:
            return ZeoppBasicPoreProperties(
                largest_included_sphere=floats[0],
                largest_free_sphere=floats[1],
                largest_included_sphere_along_free_path=floats[2],
                axis_aligned_free_sphere={"a": floats[3], "b": floats[4], "c": floats[5]},
                axis_aligned_included_sphere_along_free_path={"a": floats[6], "b": floats[7], "c": floats[8]},
            )
    raise ZeoppParseError(f"Could not parse Zeo++ -resex output: {resex_output_path}")


def _parse_chan_output(chan_output_path: Path, stdout_text: str) -> ZeoppChannelSummary:
    if not chan_output_path.is_file():
        raise ZeoppParseError(f"Zeo++ -chan output file was not created: {chan_output_path}")

    n_channels: int | None = None
    n_pockets: int | None = None
    channel_dimensionalities: tuple[int, ...] = ()
    probe_radius: float | None = None
    probe_diameter: float | None = None
    summary_lis: float | None = None
    summary_lfs: float | None = None
    summary_lisfp: float | None = None
    channel_entries: list[ZeoppChannelEntry] = []

    dimensionality_pattern = re.compile(
        r"(\d+)\s+channels identified of dimensionality\s+([0-9 ]+)\s*$",
        re.IGNORECASE,
    )
    channel_pattern = re.compile(
        r"Channel\s+(\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$",
        re.IGNORECASE,
    )
    summary_pattern = re.compile(
        r"summary\(Max_of_columns_above\)\s+"
        r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+"
        r"probe_rad:\s*(-?\d+(?:\.\d+)?)\s+probe_diam:\s*(-?\d+(?:\.\d+)?)",
        re.IGNORECASE,
    )
    identified_pattern = re.compile(r"Identified\s+(\d+)\s+channels\s+and\s+(\d+)\s+pockets\.", re.IGNORECASE)

    identified_match = identified_pattern.search(stdout_text)
    if identified_match is not None:
        n_channels = int(identified_match.group(1))
        n_pockets = int(identified_match.group(2))

    for line in chan_output_path.read_text(encoding="utf-8").splitlines():
        dimensionality_match = dimensionality_pattern.search(line)
        if dimensionality_match is not None:
            n_channels = int(dimensionality_match.group(1))
            channel_dimensionalities = tuple(int(value) for value in dimensionality_match.group(2).split())

        channel_match = channel_pattern.search(line)
        if channel_match is not None:
            channel_entries.append(
                ZeoppChannelEntry(
                    index=int(channel_match.group(1)),
                    largest_included_sphere=float(channel_match.group(2)),
                    largest_free_sphere=float(channel_match.group(3)),
                    largest_included_sphere_along_free_path=float(channel_match.group(4)),
                )
            )

        summary_match = summary_pattern.search(line)
        if summary_match is not None:
            summary_lis = float(summary_match.group(1))
            summary_lfs = float(summary_match.group(2))
            summary_lisfp = float(summary_match.group(3))
            probe_radius = float(summary_match.group(4))
            probe_diameter = float(summary_match.group(5))

    if n_channels is None:
        n_channels = len(channel_entries)
    if n_channels is None:
        raise ZeoppParseError(f"Could not parse Zeo++ -chan channel count: {chan_output_path}")
    if probe_radius is None:
        raise ZeoppParseError(f"Could not parse Zeo++ -chan probe radius: {chan_output_path}")

    if summary_lis is None and channel_entries:
        summary_lis = max(entry.largest_included_sphere for entry in channel_entries)
        summary_lfs = max(entry.largest_free_sphere for entry in channel_entries)
        summary_lisfp = max(entry.largest_included_sphere_along_free_path for entry in channel_entries)

    return ZeoppChannelSummary(
        n_channels=n_channels,
        n_pockets=n_pockets,
        channel_dimensionality=max(channel_dimensionalities) if channel_dimensionalities else None,
        channel_dimensionalities=channel_dimensionalities,
        probe_radius=probe_radius,
        probe_diameter=probe_diameter,
        largest_included_sphere=summary_lis,
        largest_free_sphere=summary_lfs,
        largest_included_sphere_along_free_path=summary_lisfp,
        channels=tuple(channel_entries),
    )


def _parse_surface_area_output(sa_output_path: Path) -> ZeoppSurfaceAreaProperties:
    if not sa_output_path.is_file():
        raise ZeoppParseError(f"Zeo++ -sa output file was not created: {sa_output_path}")

    lines = [line.rstrip() for line in sa_output_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 3:
        raise ZeoppParseError(f"Zeo++ -sa output is incomplete: {sa_output_path}")

    metrics = _parse_named_metrics(lines[0])
    channel_values = _parse_counted_values(lines[1], "Number_of_channels", "Channel_surface_area_A^2")
    pocket_values = _parse_counted_values(lines[2], "Number_of_pockets", "Pocket_surface_area_A^2")

    return ZeoppSurfaceAreaProperties(
        unitcell_volume=_require_metric(metrics, "Unitcell_volume", sa_output_path),
        density=_require_metric(metrics, "Density", sa_output_path),
        accessible_surface_area_a2=_require_metric(metrics, "ASA_A^2", sa_output_path),
        accessible_surface_area_m2_cm3=_require_metric(metrics, "ASA_m^2/cm^3", sa_output_path),
        accessible_surface_area_m2_g=_require_metric(metrics, "ASA_m^2/g", sa_output_path),
        non_accessible_surface_area_a2=_require_metric(metrics, "NASA_A^2", sa_output_path),
        non_accessible_surface_area_m2_cm3=_require_metric(metrics, "NASA_m^2/cm^3", sa_output_path),
        non_accessible_surface_area_m2_g=_require_metric(metrics, "NASA_m^2/g", sa_output_path),
        n_channels=channel_values[0],
        channel_surface_area_a2=channel_values[1],
        n_pockets=pocket_values[0],
        pocket_surface_area_a2=pocket_values[1],
    )


def _parse_volume_output(vol_output_path: Path) -> ZeoppVolumeProperties:
    if not vol_output_path.is_file():
        raise ZeoppParseError(f"Zeo++ -vol output file was not created: {vol_output_path}")

    lines = [line.rstrip() for line in vol_output_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 3:
        raise ZeoppParseError(f"Zeo++ -vol output is incomplete: {vol_output_path}")

    metrics = _parse_named_metrics(lines[0])
    channel_values = _parse_counted_values(lines[1], "Number_of_channels", "Channel_volume_A^3")
    pocket_values = _parse_counted_values(lines[2], "Number_of_pockets", "Pocket_volume_A^3")

    return ZeoppVolumeProperties(
        unitcell_volume=_require_metric(metrics, "Unitcell_volume", vol_output_path),
        density=_require_metric(metrics, "Density", vol_output_path),
        accessible_volume_a3=_require_metric(metrics, "AV_A^3", vol_output_path),
        accessible_volume_fraction=_require_metric(metrics, "AV_Volume_fraction", vol_output_path),
        accessible_volume_cm3_g=_require_metric(metrics, "AV_cm^3/g", vol_output_path),
        non_accessible_volume_a3=_require_metric(metrics, "NAV_A^3", vol_output_path),
        non_accessible_volume_fraction=_require_metric(metrics, "NAV_Volume_fraction", vol_output_path),
        non_accessible_volume_cm3_g=_require_metric(metrics, "NAV_cm^3/g", vol_output_path),
        n_channels=channel_values[0],
        channel_volume_a3=channel_values[1],
        n_pockets=pocket_values[0],
        pocket_volume_a3=pocket_values[1],
    )


def _parse_accessibility_output(axs_output_path: Path) -> ZeoppAccessibilitySummary:
    if not axs_output_path.is_file():
        raise ZeoppParseError(f"Zeo++ -axs output file was not created: {axs_output_path}")

    states = [
        line.strip().lower()
        for line in axs_output_path.read_text(encoding="utf-8").splitlines()
        if line.strip().lower() in {"true", "false"}
    ]
    if not states:
        raise ZeoppParseError(f"Could not parse Zeo++ -axs accessibility flags: {axs_output_path}")

    n_accessible = sum(state == "true" for state in states)
    n_total = len(states)
    return ZeoppAccessibilitySummary(
        n_voronoi_nodes=n_total,
        n_accessible_nodes=n_accessible,
        n_inaccessible_nodes=n_total - n_accessible,
        accessible_fraction=(n_accessible / n_total) if n_total else None,
    )


def _parse_named_metrics(text: str) -> dict[str, float]:
    return {
        match.group(1): float(match.group(2))
        for match in re.finditer(r"([A-Za-z][A-Za-z0-9_^/.]+):\s*(-?\d+(?:\.\d+)?)", text)
    }


def _require_metric(metrics: Mapping[str, float], key: str, source_path: Path) -> float:
    value = metrics.get(key)
    if value is None:
        raise ZeoppParseError(f"Could not parse {key} from {source_path}")
    return float(value)


def _parse_counted_values(text: str, count_key: str, values_key: str) -> tuple[int, tuple[float, ...]]:
    pattern = re.compile(rf"{re.escape(count_key)}:\s*(\d+)\s+{re.escape(values_key)}:\s*(.*)$")
    match = pattern.search(text)
    if match is None:
        raise ZeoppParseError(f"Could not parse {count_key} / {values_key} from Zeo++ output line: {text}")
    values = tuple(_extract_all_floats(match.group(2)))
    return int(match.group(1)), values


def _parse_trailing_floats(text: str) -> list[float]:
    match = re.search(r"(-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)+)\s*$", text)
    if match is None:
        return []
    return [float(token) for token in match.group(1).split()]


def _extract_all_floats(text: str) -> list[float]:
    return [float(match.group(0)) for match in re.finditer(r"-?\d+(?:\.\d+)?", text)]


def _format_radius_label(value: float) -> str:
    return f"{value:.3f}".replace(".", "p")


__all__ = [
    "COFKIT_ZEOPP_ENV_VAR",
    "ZeoppAccessibilitySummary",
    "ZeoppAnalysisResult",
    "ZeoppBaselineResult",
    "ZeoppBasicPoreProperties",
    "ZeoppChannelEntry",
    "ZeoppChannelSummary",
    "ZeoppConfigurationError",
    "ZeoppError",
    "ZeoppExecutionError",
    "ZeoppParseError",
    "ZeoppProbeScanResult",
    "ZeoppProbeScanSettings",
    "ZeoppSurfaceAreaProperties",
    "ZeoppVolumeProperties",
    "analyze_zeopp_basic_pore_properties",
    "analyze_zeopp_pore_properties",
    "resolve_zeopp_binary",
]
