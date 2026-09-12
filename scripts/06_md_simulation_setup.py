"""
Molecular dynamics setup and production run using OpenMM (manuscript Section
2.7), for the shortlist of top hits that survive docking + ADMET triage.

Ligand parametrization uses the OpenFF ("Sage") small-molecule force field via
`openff-toolkit` + `openmmforcefields`, avoiding any dependency on
commercial/AMBER-licensed antechamber tools -- keeps the pipeline fully
open-source, matching the paper's "open-science" framing.

Protocol (stated explicitly so it can be copied into Methods 2.7):
  1. Load receptor-ligand complex from the best-scoring Vina pose (converted
     back to a ligand SDF/mol) + prepared receptor PDB.
  2. Parametrize protein with Amber ff14SB, ligand with OpenFF Sage 2.1.
  3. Solvate in TIP3P water, 10 A padding, neutralize with Na+/Cl- to 0.15 M.
  4. Energy-minimize (L-BFGS, OpenMM default tolerance).
  5. NVT equilibration 100 ps (Langevin, 300 K), NPT equilibration 100 ps
     (Monte Carlo barostat, 1 atm).
  6. Production run (default: 100 ns; scale down and state the change if
     compute-constrained -- see critique doc Section 6, item 4).
  7. Save trajectory (DCD) + checkpoint for MM-PBSA post-processing.
"""
from __future__ import annotations

from pathlib import Path

from openff.toolkit import Molecule
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
from openmm.app import (DCDReporter, ForceField, HBonds, Modeller, PDBFile,
                         PME, Simulation, StateDataReporter)
from openmmforcefields.generators import SMIRNOFFTemplateGenerator


def build_system(receptor_pdb: str, ligand_sdf: str, out_dir: str,
                  padding_nm: float = 1.0, ionic_strength_m: float = 0.15):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    protein = PDBFile(receptor_pdb)
    ligand = Molecule.from_file(ligand_sdf)

    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand, forcefield="openff-2.1.0")
    forcefield.registerTemplateGenerator(smirnoff.generator)

    modeller = Modeller(protein.topology, protein.positions)
    modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0].to_openmm())
    modeller.addSolvent(forcefield, padding=padding_nm * unit.nanometer,
                         ionicStrength=ionic_strength_m * unit.molar,
                         neutralize=True)

    system = forcefield.createSystem(
        modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer, constraints=HBonds,
    )
    return system, modeller, out


def minimize_and_equilibrate(system, modeller, out_dir: Path,
                              temperature_k: float = 300.0,
                              nvt_steps: int = 50_000, npt_steps: int = 50_000):
    integrator = LangevinMiddleIntegrator(
        temperature_k * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds)

    platform = Platform.getPlatformByName("CUDA") if _cuda_available() else Platform.getPlatformByName("CPU")
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)
    simulation.step(nvt_steps)  # NVT equilibration

    system.addForce(MonteCarloBarostat(1 * unit.atmosphere, temperature_k * unit.kelvin))
    simulation.context.reinitialize(preserveState=True)
    simulation.step(npt_steps)  # NPT equilibration

    simulation.saveCheckpoint(str(out_dir / "equilibrated.chk"))
    return simulation


def run_production(simulation: Simulation, out_dir: Path, steps: int = 50_000_000,
                    report_interval: int = 5000):
    """Default `steps` * 2 fs = 100 ns. Reduce and document the change if
    compute-constrained (critique doc Section 6, item 4)."""
    simulation.reporters.append(DCDReporter(str(out_dir / "production.dcd"), report_interval))
    simulation.reporters.append(StateDataReporter(
        str(out_dir / "production.log"), report_interval,
        step=True, potentialEnergy=True, temperature=True, volume=True, speed=True))
    simulation.step(steps)
    simulation.saveCheckpoint(str(out_dir / "production_final.chk"))


def _cuda_available() -> bool:
    try:
        Platform.getPlatformByName("CUDA")
        return True
    except Exception:
        return False


if __name__ == "__main__":
    system, modeller, out_dir = build_system(
        receptor_pdb="data/receptors/CHK1_prepared.pdb",
        ligand_sdf="data/docking/CHK1/top_hit.sdf",
        out_dir="data/md/CHK1_top_hit",
    )
    sim = minimize_and_equilibrate(system, modeller, out_dir)
    run_production(sim, out_dir)
