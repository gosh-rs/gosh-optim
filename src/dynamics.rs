// [[file:../optim.note::899e3829][899e3829]]
use super::*;

use vecfx::*;
// 899e3829 ends here

// [[file:../optim.note::cc8bb4f6][cc8bb4f6]]
struct MoleculeDynamics<F: FnMut(&[f64], &mut [f64]) -> Result<f64>> {
    dynamics: crate::potential::Dynamics<F>,

    mass: Vec<f64>,
    velocity: Vec<f64>,
}
// cc8bb4f6 ends here

// [[file:../optim.note::e0fc00ce][e0fc00ce]]
impl<F> MoleculeDynamics<F>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    /// update velocity and positions in Velocity Verlet Algorithm
    fn velocity_verlet_update(&mut self, dt: f64) -> Result<()> {
        let m = self.mass.as_vector_slice();
        let v = self.velocity.as_vector_slice();
        let r = self.dynamics.positions().as_vector_slice();
        let f = self.dynamics.get_forces()?.as_vector_slice();

        // update positions
        let dr = v * dt + 0.5 * f.component_div(&m) * dt.powi(2);
        self.dynamics.step_toward(dr.as_slice());

        // update velecities
        let f_new = self.dynamics.get_forces()?.as_vector_slice();
        let v_new = dr / dt + 0.5 * f_new.component_div(&m) * dt;
        self.velocity.copy_from_slice(v_new.as_slice());

        Ok(())
    }
}
// e0fc00ce ends here

// [[file:../optim.note::0f175760][0f175760]]
impl<F> MoleculeDynamics<F>
where
    F: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    /// Trajectory Propagation
    pub fn propagate(&mut self, timestep: f64) -> Result<()> {
        self.velocity_verlet_update(timestep)
    }
}
// 0f175760 ends here
