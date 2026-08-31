"""Diagnostic-only Phase-18B export of the pre-modal UVLM contract."""

import hashlib
import json
import os

import h5py
import numpy as np
import scipy.sparse as sp

from sharpy.utils.solver_interface import BaseSolver, solver
import sharpy.utils.settings as settings_utils


def _dense(value):
    if sp.issparse(value):
        return value.toarray()
    return np.asarray(value)


def _write_state_space(group, state_space):
    group.create_dataset("A", data=_dense(state_space.A))
    group.create_dataset("B", data=_dense(state_space.B))
    group.create_dataset("C", data=_dense(state_space.C))
    group.create_dataset("D", data=_dense(state_space.D))
    group.attrs["dt"] = np.nan if state_space.dt is None else float(state_space.dt)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@solver
class Phase18PremodalExtractor(BaseSolver):
    """Save matrices immediately surrounding structural modal projection."""

    solver_id = "Phase18PremodalExtractor"
    solver_classification = "post-processor"

    settings_types = {"enabled": "bool", "subfolder": "str"}
    settings_default = {"enabled": False, "subfolder": "phase18_premodal"}
    settings_description = {
        "enabled": "Enable the diagnostic export. Off by default.",
        "subfolder": "Output subfolder below the SHARPy case output.",
    }

    def __init__(self):
        self.data = None
        self.settings = None

    def initialise(self, data, custom_settings=None, restart=False):
        del restart
        self.data = data
        self.settings = (data.settings[self.solver_id]
                         if custom_settings is None else custom_settings)
        settings_utils.to_custom_types(
            self.settings, self.settings_types, self.settings_default)

    def run(self, **kwargs):
        del kwargs
        if not self.settings["enabled"]:
            return self.data

        linear = self.data.linear.linear_system
        uvlm = linear.uvlm
        beam = linear.beam
        krylov = uvlm.rom["Krylov"]

        required_couplings = ("Kas", "Ksa", "in_mode_gain", "out_mode_gain")
        missing = [name for name in required_couplings
                   if name not in linear.couplings]
        if missing:
            raise RuntimeError(
                "Phase-18 pre-modal extraction missing couplings: "
                + ", ".join(missing))

        folder = os.path.join(self.data.output_folder,
                              str(self.settings["subfolder"]))
        os.makedirs(folder, exist_ok=True)
        artifact = os.path.join(folder, "premodal_contract.h5")

        with h5py.File(artifact, "w") as handle:
            handle.attrs["schema_version"] = 2
            handle.attrs["case"] = self.data.settings["SHARPy"]["case"]
            handle.attrs["projection_order"] = (
                "physical_discrete->nodal_Kas_Ksa->modal->Krylov")
            handle.attrs["enabled_by_explicit_option"] = True

            _write_state_space(handle.create_group("physical_discrete"),
                               uvlm.sys.SS)
            augmented = handle.create_group("physical_augmented_discrete")
            augmented.create_dataset("A", data=_dense(krylov.ss.A))
            augmented.create_dataset("B", data=_dense(uvlm.B_to_vertex_forces))
            augmented.create_dataset("C", data=_dense(uvlm.C_to_vertex_forces))
            augmented.create_dataset("D", data=_dense(uvlm.D_to_vertex_forces))
            augmented.attrs["dt"] = float(krylov.ss.dt)
            _write_state_space(handle.create_group("modal_pre_krylov"),
                               krylov.ss)
            _write_state_space(handle.create_group("krylov_discrete"),
                               krylov.ssrom)
            _write_state_space(handle.create_group("legacy_saved_uvlm"),
                               uvlm.ss)

            maps = handle.create_group("maps")
            maps.create_dataset("Kas_nodal_to_aero",
                                data=_dense(linear.couplings["Kas"].value))
            maps.create_dataset("Ksa_aero_to_nodal",
                                data=_dense(linear.couplings["Ksa"].value))
            maps.create_dataset("modal_input",
                                data=_dense(linear.couplings["in_mode_gain"].value))
            maps.create_dataset("modal_output",
                                data=_dense(linear.couplings["out_mode_gain"].value))
            maps.create_dataset("structural_modes", data=_dense(beam.sys.U))
            maps.create_dataset("Kdisp", data=_dense(linear.Kdisp))
            maps.create_dataset("Kdisp_vel", data=_dense(linear.Kdisp_vel))
            maps.create_dataset("Kvel_disp", data=_dense(linear.Kvel_disp))
            maps.create_dataset("Kvel_vel", data=_dense(linear.Kvel_vel))
            maps.create_dataset("Kforces", data=_dense(linear.Kforces))

            # SHARPy assembles Kforces with the flexible nodal rows first,
            # followed by the rigid resultant force and moment rows.  Ksa
            # deliberately retains only the flexible rows for structural
            # propagation; retaining the next six rows preserves the direct
            # aerodynamic resultant without reconstructing it from modal
            # generalized forces.
            kforces = _dense(linear.Kforces)
            n_flexible = int(beam.sys.num_dof)
            if kforces.shape[0] < n_flexible + 6:
                raise RuntimeError(
                    "Phase-18 root output requires six rigid-resultant "
                    "rows after the flexible structural force rows")
            root_map = kforces[n_flexible:n_flexible + 6, :]
            maps.create_dataset("Kforces_full", data=kforces)
            maps.create_dataset("root_resultant_from_vertex", data=root_map)

            physical_output = handle.create_group("physical_output")
            nodal_c = _dense(linear.couplings["Ksa"].value).dot(
                _dense(uvlm.C_to_vertex_forces))
            nodal_d = _dense(linear.couplings["Ksa"].value).dot(
                _dense(uvlm.D_to_vertex_forces))
            root_c = root_map.dot(_dense(uvlm.C_to_vertex_forces))
            root_d = root_map.dot(_dense(uvlm.D_to_vertex_forces))
            physical_output.create_dataset("nodal_C", data=nodal_c)
            physical_output.create_dataset("nodal_D", data=nodal_d)
            physical_output.create_dataset("root_C", data=root_c)
            physical_output.create_dataset("root_D", data=root_d)
            physical_output.create_dataset(
                "steady_nodal_force",
                data=_dense(linear.couplings["Ksa"].value).dot(
                    np.asarray(uvlm.linearisation_vectors["forces_aero"])))
            physical_output.create_dataset(
                "steady_root_resultant",
                data=root_map.dot(np.asarray(
                    uvlm.linearisation_vectors["forces_aero"])))
            physical_output.attrs["frame"] = "SHARPy structural A frame"
            physical_output.attrs["reference_point"] = "structural root"
            physical_output.attrs["units"] = "active SHARPy force scaling"
            physical_output.attrs["component_order"] = "Fx,Fy,Fz,Mx,My,Mz"
            physical_output.attrs["multiplicity"] = "full generated wing model"
            physical_output.attrs["sign_convention"] = (
                "positive applied aerodynamic resultant on structure")
            physical_output.attrs["ownership"] = (
                "external aerodynamic load only; excludes structural "
                "internal, preload, inertia, gravity, thrust and tail loads")

            projectors = handle.create_group("krylov_projectors")
            projectors.create_dataset("V", data=_dense(krylov.V.value
                                                        if hasattr(krylov.V, "value")
                                                        else krylov.V))
            projectors.create_dataset("W", data=_dense(krylov.W.value
                                                        if hasattr(krylov.W, "value")
                                                        else krylov.W))

            scaling = handle.create_group("scaling")
            scaling.attrs["uvlm_scaled"] = bool(uvlm.scaled)
            for name, value in uvlm.sys.ScalingFacts.items():
                scaling.attrs[name] = float(value)

            geometry = handle.create_group("geometry")
            for name in ("zeta", "zeta_dot", "u_ext", "forces_aero"):
                geometry.create_dataset(name, data=np.asarray(
                    uvlm.linearisation_vectors[name]))

            static = handle.create_group("static_equilibrium")
            undeformed = self.data.structure.ini_info
            equilibrium = self.data.structure.timestep_info[0]
            static.create_dataset("undeformed_nodal_position",
                                  data=np.asarray(undeformed.pos))
            static.create_dataset("deformed_nodal_position",
                                  data=np.asarray(equilibrium.pos))
            static.create_dataset("undeformed_element_crv",
                                  data=np.asarray(undeformed.psi))
            static.create_dataset("deformed_element_crv",
                                  data=np.asarray(equilibrium.psi))
            static.create_dataset("frame_quaternion",
                                  data=np.asarray(equilibrium.quat))
            static.create_dataset("structural_q",
                                  data=np.asarray(equilibrium.q))
            static.create_dataset("aerodynamic_nodal_force",
                                  data=np.asarray(equilibrium.postproc_node[
                                      "aero_steady_forces"]))
            static.create_dataset("gravity_nodal_force",
                                  data=np.asarray(equilibrium.gravity_forces))
            static.create_dataset("constraint_nodal_reaction",
                                  data=np.asarray(
                                      equilibrium.forces_constraints_nodes))
            if "strain" not in equilibrium.postproc_cell or \
                    "loads" not in equilibrium.postproc_cell:
                raise RuntimeError(
                    "Phase-18 static export requires BeamLoads before the "
                    "dynamic step")
            static.create_dataset("sectional_strain_curvature",
                                  data=np.asarray(
                                      equilibrium.postproc_cell["strain"]))
            static.create_dataset("sectional_internal_resultant",
                                  data=np.asarray(
                                      equilibrium.postproc_cell["loads"]))
            control = np.asarray(
                self.data.aero.data_dict["control_surface_deflection"])
            static.create_dataset("control_surface_deflection", data=control)
            external = (np.asarray(equilibrium.postproc_node[
                "aero_steady_forces"]) +
                np.asarray(equilibrium.gravity_forces))
            root = np.asarray(undeformed.pos)[np.flatnonzero(
                np.asarray(self.data.structure.boundary_conditions) == 1)[0]]
            support = -np.sum(external, axis=0)
            support[3:6] = -np.sum(
                external[:, 3:6] +
                np.cross(np.asarray(equilibrium.pos) - root,
                         external[:, 0:3]), axis=0)
            static.create_dataset("support_reaction", data=support)
            static.attrs["frame"] = "SHARPy structural A frame"
            static.attrs["root_reference"] = "clamped structural node"
            static.attrs["multiplicity"] = "full generated wing; no mirroring"
            static.attrs["force_units"] = "N"
            static.attrs["moment_units"] = "N*m"
            static.attrs["support_reaction_derivation"] = (
                "negative static aerodynamic-plus-gravity nodal wrench "
                "shifted to the clamped root")

        manifest = {
            "schema_version": 2,
            "artifact": os.path.basename(artifact),
            "sha256": _sha256(artifact),
            "case": self.data.settings["SHARPy"]["case"],
            "source_dt": float(uvlm.sys.SS.dt),
            "dimensions": {
                "physical_states": int(uvlm.sys.SS.states),
                "physical_inputs": int(uvlm.sys.SS.inputs),
                "physical_outputs": int(uvlm.sys.SS.outputs),
                "augmented_physical_states": int(krylov.ss.states),
                "augmented_physical_inputs": int(uvlm.B_to_vertex_forces.shape[1]),
                "augmented_physical_outputs": int(uvlm.C_to_vertex_forces.shape[0]),
                "modal_inputs": int(krylov.ss.inputs),
                "modal_outputs": int(krylov.ss.outputs),
                "krylov_states": int(krylov.ssrom.states),
            },
        }
        manifest_path = os.path.join(folder, "premodal_manifest.json")
        with open(manifest_path, "w", encoding="utf-8") as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")

        print(f"[P18-PREMODAL] wrote {artifact}")
        return self.data
