// [[file:../optim.note::a197ff17][a197ff17]]
use super::*;
use crate::vars::Vars;

use fire::fire;
use lbfgs::lbfgs_iter;
// a197ff17 ends here

// [[file:../optim.note::*base][base:1]]
#[derive(Debug, Clone)]
/// A helper struct containing information on optimization step.
pub struct OptimizedIter {
    /// The number of calls for potential evaluation.
    pub ncalls: usize,
    /// Current fmax criterion of forces in optimization.
    pub fmax: f64,
    /// Current energy in optimization.
    pub energy: f64,
}
// base:1 ends here

// [[file:../optim.note::fe25e584][fe25e584]]
pub fn optimize_geometry_iter<'a>(potential: &'a mut Dynamics) -> Box<dyn Iterator<Item = OptimizedIter> + 'a> {
    let vars = Vars::from_env();
    if vars.algorithm == "FIRE" {
        info!("Optimizing using FIRE algorithm ...");
        let x_init = potential.position().to_vec();
        let mut opt = fire()
            .with_max_step(vars.max_step_size)
            .with_max_cycles(vars.max_evaluations);
        let steps = opt.minimize_iter(x_init, move |x: &[f64], o: &mut fire::Output| {
            potential.set_position(x);
            let energy = potential.get_energy()?;
            let force = potential.get_force()?;
            o.fx = energy;
            o.gx.vecncpy(force);
            let fmax = force.iter().map(|x| x.abs()).float_max();
            let ncalls = potential.ncalls();
            let progress = OptimizedIter { ncalls, energy, fmax };
            Ok(progress)
        });
        Box::new(steps.map(|progress| progress.extra))
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

        let x_init = potential.position().to_vec();
        let steps = opt.minimize(x_init, move |x: &[f64], o: &mut lbfgs::Output| {
            potential.set_position(x);
            let energy = potential.get_energy()?;
            let force = potential.get_force()?;
            o.fx = energy;
            o.gx.vecncpy(force);
            let fmax = force.iter().map(|x| x.abs()).float_max();
            let ncalls = potential.ncalls();
            let progress = OptimizedIter { ncalls, energy, fmax };
            Ok(progress)
        }).expect("optimize lbfgs");
        Box::new(steps.map(|progress| progress.extra))
    }
}
// fe25e584 ends here