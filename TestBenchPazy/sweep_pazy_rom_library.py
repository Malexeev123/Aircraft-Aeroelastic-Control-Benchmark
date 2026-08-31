"""Generate Pazy SHARPy source points for the MATLAB ROM workflow.

The script calls the installed SHARPy/XBeam interfaces without modifying
either dependency. It defaults to a plan preview; pass ``--execute`` to run.
Each point is preserved under ``library_source`` unless ``--overwrite`` is
stated explicitly.
"""

from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import scipy.io

from get_settings_krylov import get_settings_krylov
from get_settings_udp import get_settings_udp


CASE_NAME_ROM = "pazy_krylov_ROM"
CASE_NAME_OPEN_LOOP = "pazy_open_loop_gust_response"


def load_sharpy():
    try:
        import sharpy.cases.templates.flying_wings as flying_wings
        import sharpy.sharpy_main as sharpy_main
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "SHARPy is not available in this Python environment. Activate "
            "the established WSL SHARPy environment before using --execute."
        ) from error
    return flying_wings, sharpy_main


@dataclass(frozen=True)
class GenerationSettings:
    density: float = 1.225
    chordwise_panels: int = 4
    spanwise_nodes: int = 16
    wake_length_factor: int = 10
    control_surfaces: int = 2
    modes: int = 10
    cores: int = 4
    simulation_time: float = 1.0
    gust_length: float = 10.0
    gust_intensity: float = 0.02

    @property
    def gust(self) -> dict[str, float | str]:
        return {
            "gust_shape": "1-cos",
            "gust_length": self.gust_length,
            "gust_intensity": self.gust_intensity,
            "gust_offset": 0.0,
        }


def point_name(speed: float, alpha_deg: float) -> str:
    speed_token = f"{int(round(speed)):03d}"
    sign = "m" if alpha_deg < 0 else "p"
    alpha_token = f"{abs(int(round(alpha_deg))):02d}"
    return f"pt_U{speed_token}_alpha_{sign}{alpha_token}"


def parse_values(values: Iterable[float] | None,
                 default: list[float]) -> list[float]:
    parsed = default if values is None else [float(value) for value in values]
    if not parsed or not all(np.isfinite(parsed)):
        raise ValueError("Every grid coordinate must be finite.")
    return parsed


def make_model(root: Path, case_name: str, speed: float, alpha_deg: float,
               settings: GenerationSettings, physical_time: float):
    flying_wings, _ = load_sharpy()
    return flying_wings.PazyControlSurface(
        M=settings.chordwise_panels,
        N=settings.spanwise_nodes,
        Mstar_fact=settings.wake_length_factor,
        u_inf=speed,
        alpha=alpha_deg,
        rho=settings.density,
        n_surfaces=2,
        route=str(root / "cases" / case_name),
        case_name=case_name,
        physical_time=physical_time,
    )


def prepare_model_files(model, fixed_dt: float | None = None) -> None:
    model.clean_test_files()
    if fixed_dt is not None:
        model.dt_factor = fixed_dt * model.M * model.u_inf / model.main_chord
        model.physical_time = fixed_dt
    model.update_derived_params()
    if fixed_dt is not None and abs(model.dt - fixed_dt) > 1e-13:
        raise RuntimeError(
            f"SHARPy resolved dt={model.dt:.16e}; expected {fixed_dt:.16e}."
        )
    model.generate_aero_file()
    model.generate_fem_file()


def run_open_loop_reference(root: Path, speed: float, alpha_deg: float,
                            settings: GenerationSettings) -> np.ndarray:
    _, sharpy_main = load_sharpy()
    model = make_model(root, CASE_NAME_OPEN_LOOP, speed, alpha_deg, settings,
                       settings.simulation_time)
    prepare_model_files(model)
    flow = [
        "BeamLoader", "AerogridLoader", "StaticCoupled",
        "WriteVariablesTime", "DynamicCoupled", "AerogridPlot", "BeamPlot",
    ]
    model.set_default_config_dict()
    model.config = get_settings_udp(
        model, flow, num_cores=settings.cores,
        wake_length_factor=settings.wake_length_factor,
        output_folder=str(root / "output"), gust=True,
        gust_settings=settings.gust,
    )
    model.config.write()
    config_path = Path(model.route) / f"{model.case_name}.sharpy"
    sharpy_main.main(["", str(config_path)])

    tip_file = (
        root / "output" / model.case_name / "WriteVariablesTime" /
        f"struct_pos_node{model.num_node_surf}.dat"
    )
    if not tip_file.is_file():
        raise FileNotFoundError(f"SHARPy did not write the tip history: {tip_file}")
    tip_z = np.loadtxt(tip_file, ndmin=2)[:, -1]
    return 100.0 * tip_z / (0.5 * model.b_ref)


def run_rom_source(root: Path, speed: float, alpha_deg: float,
                   settings: GenerationSettings, fixed_dt: float | None,
                   extract_premodal: bool):
    _, sharpy_main = load_sharpy()
    probe = make_model(root, CASE_NAME_ROM, speed, alpha_deg, settings, 1.0)
    probe.update_derived_params()
    point_dt = probe.dt if fixed_dt is None else fixed_dt
    model = make_model(root, CASE_NAME_ROM, speed, alpha_deg, settings, point_dt)
    prepare_model_files(model, fixed_dt=fixed_dt)

    flow = ["BeamLoader", "AerogridLoader", "StaticCoupled"]
    if extract_premodal:
        # Importing registers this project-owned postprocessor with SHARPy.
        import phase18_premodal_extractor  # noqa: F401

        flow.append("BeamLoads")
    flow.extend(["WriteVariablesTime", "DynamicCoupled", "Modal", "LinearAssembler"])
    if extract_premodal:
        flow.append("Phase18PremodalExtractor")
    flow.extend(["SaveData", "SaveParametricCase"])

    rom_settings = {
        "use": True,
        "rom_method": "Krylov",
        "rom_method_settings": {
            "Krylov": {
                "algorithm": "mimo_rational_arnoldi",
                "r": 4,
                "frequency": np.array([0]),
                "single_side": "observability",
            }
        },
    }
    model.set_default_config_dict()
    model.config = get_settings_krylov(
        model, flow, num_cores=settings.cores,
        wake_length_factor=settings.wake_length_factor,
        output_folder=str(root / "output"), gust=True,
        gust_settings=settings.gust, num_modes=settings.modes,
        rom_settings=rom_settings, remove_gust_input_in_statespace=False,
        unsteady_force_distribution=False, gravity_on=True,
    )
    if extract_premodal:
        model.config["Phase18PremodalExtractor"] = {
            "enabled": True,
            "subfolder": "phase18_premodal",
        }
    model.config.write()

    config_path = Path(model.route) / f"{model.case_name}.sharpy"
    print(
        f"[SHARPy] {point_name(speed, alpha_deg)}: "
        f"U={speed:g} m/s, alpha={alpha_deg:g} deg, dt={model.dt:.9e} s"
    )
    sharpy_main.main(["", str(config_path)])
    return model


def write_matlab_parameters(root: Path, model, speed: float, alpha_deg: float,
                            settings: GenerationSettings,
                            tip_reference: np.ndarray | None) -> None:
    destination = root / "output" / model.case_name / "linear_results"
    destination.mkdir(parents=True, exist_ok=True)
    sample_count = int(round(settings.simulation_time / model.dt)) + 1
    if tip_reference is None:
        tip_reference = np.zeros(sample_count)
    tip_reference = np.asarray(tip_reference, dtype=float).reshape(-1)
    if tip_reference.size < sample_count:
        pad = tip_reference[-1] if tip_reference.size else 0.0
        tip_reference = np.pad(
            tip_reference, (0, sample_count - tip_reference.size),
            constant_values=pad,
        )
    tip_reference = tip_reference[:sample_count]

    scipy.io.savemat(
        destination / "simulation_parameters.mat",
        {
            "num_aero_states": 4 * settings.modes,
            "num_modes": settings.modes,
            "u_inf": float(speed),
            "simulation_time": settings.simulation_time,
            "gust_length": settings.gust_length,
            "gust_intensity": settings.gust_intensity,
            "num_control_surfaces": settings.control_surfaces,
            "n_nodes": settings.spanwise_nodes,
            "control_input_start": 194,
            "gust_input_start": 193,
            "aoa_deg": float(alpha_deg),
            "normalised_tip_displacement": tip_reference,
            "dt": float(model.dt),
            "dt_factor": float(model.dt_factor),
        },
        do_compression=True,
    )


def snapshot_point(root: Path, point_id: str, overwrite: bool) -> Path:
    destination = root / "library_source" / CASE_NAME_ROM / point_id
    if destination.exists():
        if not overwrite:
            raise FileExistsError(
                f"Source point already exists: {destination}. "
                "Use --overwrite only when replacement is intentional."
            )
        shutil.rmtree(destination)
    (destination / "cases").mkdir(parents=True)
    (destination / "output").mkdir(parents=True)
    shutil.copytree(root / "cases" / CASE_NAME_ROM,
                    destination / "cases" / CASE_NAME_ROM)
    shutil.copytree(root / "output" / CASE_NAME_ROM,
                    destination / "output" / CASE_NAME_ROM)
    return destination


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--speed", type=float, action="append",
                        help="Airspeed in m/s; repeat for multiple points.")
    parser.add_argument("--alpha", type=float, action="append",
                        help="Angle of attack in degrees; repeat as needed.")
    parser.add_argument("--execute", action="store_true",
                        help="Run SHARPy. Otherwise print the plan only.")
    parser.add_argument("--open-loop-reference", action="store_true",
                        help="Generate the matched nonlinear wingtip reference.")
    parser.add_argument("--extract-premodal", action="store_true",
                        help="Export the project-owned premodal contract.")
    parser.add_argument("--fixed-dt", type=float,
                        help="Optional common positive source timestep in seconds.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace an existing point snapshot.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent / "generated",
        help=(
            "Generated-data workspace (default: TestBenchPazy/generated). "
            "Keep this separate from the tracked runtime payload."
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    speeds = parse_values(args.speed, [40.0])
    alphas = parse_values(args.alpha, [1.0])
    if any(speed <= 0 for speed in speeds):
        raise ValueError("Airspeed must be positive.")
    if args.fixed_dt is not None and args.fixed_dt <= 0:
        raise ValueError("--fixed-dt must be positive.")

    root = args.root.resolve()
    points = [(speed, alpha) for speed in speeds for alpha in alphas]
    print(f"Pazy source-generation root: {root}")
    for speed, alpha in points:
        print(f"  {point_name(speed, alpha):24s}  U={speed:g} m/s  alpha={alpha:g} deg")
    if not args.execute:
        print("Plan only. Re-run with --execute to call SHARPy.")
        return 0

    root.mkdir(parents=True, exist_ok=True)
    settings = GenerationSettings()
    for speed, alpha in points:
        reference = None
        if args.open_loop_reference:
            reference = run_open_loop_reference(root, speed, alpha, settings)
        model = run_rom_source(
            root, speed, alpha, settings, args.fixed_dt, args.extract_premodal,
        )
        write_matlab_parameters(root, model, speed, alpha, settings, reference)
        destination = snapshot_point(
            root, point_name(speed, alpha), args.overwrite,
        )
        print(f"[complete] {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
