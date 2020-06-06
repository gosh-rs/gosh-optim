// [[file:~/Workspace/Programming/gosh-rs/optim/optim.note::*imports][imports:1]]
use gosh_core::*;

use gosh_model::{ChemicalModel, ModelProperties};

use gchemol::Molecule;
use gut::prelude::*;
use vecfx::*;
// imports:1 ends here

// [[file:~/Workspace/Programming/gosh-rs/optim/optim.note::*base][base:1]]
/// A generic interface for geometry optimization of Molecule.
pub struct Optimizer {
    fmax: f64,
    nmax: usize,
}

impl Default for Optimizer {
    fn default() -> Self {
        Self { fmax: 0.1, nmax: 100 }
    }
}

/// A helper struct containing information on optimization.
pub struct Optimized {
    /// The number of iterations in Optimzation loop.
    pub niter: usize,
    /// Final computed properties in ChemicalModel
    pub computed: ModelProperties,
}
// base:1 ends here

// [[file:~/Workspace/Programming/gosh-rs/optim/optim.note::*pub][pub:1]]
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
    pub fn optimize_geometry<M: ChemicalModel>(&self, mol: &mut Molecule, model: &mut M) -> Result<Optimized> {
        let coords = mol.positions().collect_vec().concat();
        let mask = mol.freezing_coords_mask();
        let mut x_init_masked = mask.apply(&coords);
        let mut computed = None;
        let mut opt = lbfgs::lbfgs()
            .with_max_evaluations(100)
            .with_initial_step_size(1.0 / 75.0)
            .with_max_step_size(0.1)
            .with_gradient_only()
            .with_damping(true)
            .with_max_linesearch(1)
            .with_linesearch_gtol(0.999)
            .build(&mut x_init_masked, |x_masked, gx_masked| {
                let positions = mask.unmask(x_masked, 0.0).as_3d().to_owned();
                mol.set_positions(positions);
                let mp = model.compute(&mol)?;
                let energy = mp.get_energy().ok_or(format_err!("opt: no energy"))?;
                let forces = mp.get_forces().ok_or(format_err!("opt: no forces"))?;
                let forces = mask.apply(forces.as_flat());
                gx_masked.vecncpy(&forces);

                // save for returning
                computed = Some(mp);
                Ok(energy)
            })?;

        let mut niter = 0;
        for i in 0..self.nmax {
            let state = opt.propagate()?;
            let fmax = state.gx.chunks(3).map(|v| v.vec2norm()).float_max();
            println!("iter {:4}\tEnergy = {:-12.3}\tfmax={}", i, state.fx, fmax);
            if fmax < self.fmax {
                niter = i;
                break;
            }
        }

        let mp = computed.ok_or(format_err!("model was not computed"))?;
        let optimized = Optimized { niter, computed: mp };

        Ok(optimized)
    }
}
// pub:1 ends here

// [[file:~/Workspace/Programming/gosh-rs/optim/optim.note::*test][test:1]]
#[test]
fn test_opt() -> Result<()> {
    use gchemol::prelude::*;
    use gosh_model::LennardJones;

    let filename = "tests/files/LennardJones/LJ38r.xyz";
    let mut mol = Molecule::from_file(filename)?;
    let mut lj = LennardJones::default();
    lj.derivative_order = 1;

    let _ = Optimizer::default().optimize_geometry(&mut mol, &mut lj)?;

    Ok(())
}
// test:1 ends here
