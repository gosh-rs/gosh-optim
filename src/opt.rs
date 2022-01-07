// [[file:../optim.note::*imports][imports:1]]
use gosh_core::*;

use gosh_model::{ChemicalModel, ModelProperties};

use gchemol::Molecule;
use gut::prelude::*;
use vecfx::*;
// imports:1 ends here

// [[file:../optim.note::*base][base:1]]
use gosh_database::CheckpointDb;

/// A generic interface for geometry optimization of Molecule.
pub struct Optimizer {
    fmax: f64,
    nmax: usize,
    ckpt: Option<CheckpointDb>,
    vars: crate::vars::Vars,
}

impl Default for Optimizer {
    fn default() -> Self {
        Self {
            fmax: 0.1,
            nmax: 100,
            ckpt: None,
            vars: crate::vars::Vars::from_env(),
        }
    }
}

impl Optimizer {
    /// New optimizer with max force (fmax) and max step (nmax) in iterations.
    pub fn new(fmax: f64, nmax: usize) -> Self {
        Self {
            fmax,
            nmax,
            ..Self::default()
        }
    }

    /// Set checkpoint for resuming optimization later
    pub fn checkpoint(mut self, ckpt: CheckpointDb) -> Self {
        self.ckpt = ckpt.into();
        self
    }
}

/// A helper struct containing information on optimization.
pub struct Optimized {
    /// The number of iterations in optimzation loop.
    pub niter: usize,
    /// Final fmax criterion for forces.
    pub fmax: f64,
    /// Final computed properties in ChemicalModel.
    pub computed: ModelProperties,
}
// base:1 ends here

// [[file:../optim.note::*pub/trait][pub/trait:1]]
/// A helper struct represents the output data required for molecular geometry
/// optimization.
pub struct Output {
    pub energy: Option<f64>,
    pub forces: Option<Vec<[f64; 3]>>,
}

pub trait OptimizeMolecule<U> {
    fn evaluate(&mut self, mol: &Molecule, out: &mut Output) -> Result<U>;
}

impl<T> OptimizeMolecule<ModelProperties> for T
where
    T: ChemicalModel,
{
    fn evaluate(&mut self, mol: &Molecule, out: &mut Output) -> Result<ModelProperties> {
        trace!("opt: evaluate PES");
        let mut mp = self.compute(&mol)?;

        // save for returning
        // make sure `ModelProperties` contains correct version of `Molecule`
        mp.set_molecule(mol.clone());
        out.energy = mp.get_energy();
        out.forces = mp.get_forces().cloned();

        Ok(mp)
    }
}
// pub/trait:1 ends here

// [[file:../optim.note::b17504d6][b17504d6]]
#[derive(Debug, Clone)]
/// A helper struct containing information on optimization step.
pub struct OptimizedIter<U> {
    /// The number of calls for potential evaluation.
    pub ncalls: usize,
    /// Current fmax criterion of forces in optimization.
    pub fmax: f64,
    /// Current energy in optimization.
    pub energy: f64,
    /// Extra data returned from user defined OptimizeMolecule trait method
    pub extra: U,
}

/// Optimize geometry of `mol` in potential provided by `model` (iterator version).
///
/// # Parameters
///
/// * mol: target molecule
/// * model: chemical model for evaluation energy and forces of `mol` at new positions.
///
/// # Return
///
/// Returns an iterator over optimization steps
pub fn optimize_geometry_iter<'a, M, U: 'a>(
    mol: &'a mut Molecule,
    model: &'a mut M,
) -> Box<dyn Iterator<Item = OptimizedIter<U>> + 'a>
where
    M: OptimizeMolecule<U>,
{
    let vars = crate::vars::Vars::from_env();
    dbg!(&vars);
    let coords = mol.positions().collect_vec().concat();
    let mask = mol.freezing_coords_mask();
    let mut x_init_masked = mask.apply(&coords);

    if vars.algorithm == "FIRE" {
        info!("Optimizing using FIRE algorithm ...");
        let mut opt = fire::fire()
            .with_max_step(vars.max_step_size)
            .with_max_cycles(vars.max_evaluations);

        let steps = opt.minimize_iter(x_init_masked, move |x_masked: &[f64], o_masked: &mut fire::Output| {
            let positions = mask.unmask(x_masked, 0.0).as_3d().to_owned();
            mol.update_positions(positions);
            let mut out = Output {
                energy: None,
                forces: None,
            };
            let extra = model.evaluate(&mol, &mut out)?;
            let energy = out.energy.expect("evaluate: forget to set energy?");
            let forces = out.forces.as_ref().expect("evaluate: forget to set forces?");
            let forces = mask.apply(forces.as_flat());
            trace!("opt: evaluate PES");

            o_masked.gx.vecncpy(&forces);
            o_masked.fx = energy;

            let fmax = forces.chunks(3).map(|v| v.vec2norm()).float_max();
            Ok((fmax, extra))
        });

        Box::new(steps.map(|progress| {
            let (fmax, extra) = progress.extra;
            OptimizedIter {
                fmax,
                extra,
                ncalls: progress.ncalls,
                energy: progress.fx,
            }
        }))
    } else {
        info!("Optimizing using L-BFGS algorithm ...");
        let mut opt = lbfgs::lbfgs_iter()
            .with_max_evaluations(vars.max_evaluations)
            .with_initial_step_size(vars.initial_step_size)
            .with_max_step_size(vars.max_step_size)
            .with_max_linesearch(vars.max_linesearch)
            .with_gradient_only()
            .with_damping(true)
            .with_linesearch_gtol(0.999);

        let steps = opt
            .minimize(x_init_masked, move |x_masked: &[f64], o_masked: &mut lbfgs::Output| {
                let positions = mask.unmask(x_masked, 0.0).as_3d().to_owned();
                mol.update_positions(positions);
                let mut out = Output {
                    energy: None,
                    forces: None,
                };
                let extra = model.evaluate(&mol, &mut out)?;
                let energy = out.energy.expect("evaluate: forget to set energy?");
                let forces = out.forces.as_ref().expect("evaluate: forget to set forces?");
                let forces = mask.apply(forces.as_flat());
                trace!("opt: evaluate PES");

                o_masked.gx.vecncpy(&forces);
                o_masked.fx = energy;

                let fmax = forces.chunks(3).map(|v| v.vec2norm()).float_max();
                Ok((fmax, extra))
            })
            .expect("optimize_geometry_iter");

        Box::new(steps.map(|progress| {
            let (fmax, extra) = progress.extra;
            OptimizedIter {
                fmax,
                extra,
                ncalls: progress.ncalls,
                energy: progress.fx,
            }
        }))
    }
}
// b17504d6 ends here

// [[file:../optim.note::315bd793][315bd793]]
impl Optimizer {
    /// Optimize geometry of `mol` in potential provided by `model`.
    ///
    /// # Parameters
    ///
    /// * mol: target molecule
    /// * model: chemical model for evaluation energy and forces of `mol` at new positions.
    ///
    /// # Return
    ///
    /// Returns the computed `ModelProperties` on success in final geometry.
    ///
    pub fn optimize_geometry<M: ChemicalModel>(&self, mol: &mut Molecule, model: &mut M) -> Result<Optimized> {
        // restore Molecule from ckpt
        if let Some(ckpt) = &self.ckpt {
            ckpt.restore(mol).context("restore optimized molecule from ckpt")?;
        }

        let steps = self::optimize_geometry_iter(mol, model);

        let mut computed = None;
        let mut niter = 0;
        let mut fmax = std::f64::NAN;
        for (progress, i) in steps.take(self.nmax).zip(1..) {
            // checkpointing
            if let Some(ckpt) = &self.ckpt {
                let mol = progress.extra.get_molecule().expect("no mol in mp");
                ckpt.commit(mol);
            }

            niter = i;
            fmax = progress.fmax;
            computed = progress.extra.into();
            println!("iter {:4}\tEnergy = {:-12.4}\tfmax={}", i, progress.energy, fmax);
            if fmax < self.fmax {
                info!("forces converged: {}", fmax);
                break;
            }
        }

        // FIXME: it is better to use `OptimizedIter`?
        let mp = computed.ok_or(format_err!("model was not computed"))?;
        let optimized = Optimized {
            niter,
            fmax,
            computed: mp,
        };

        Ok(optimized)
    }
}
// 315bd793 ends here
