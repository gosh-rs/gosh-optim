// [[file:../optim.note::899e3829][899e3829]]
use super::*;
// 899e3829 ends here

// [[file:../optim.note::e0fc00ce][e0fc00ce]]
struct MoleculeDynamics {
    forces: Vec<f64>,
    masses: Vec<f64>,
    positions: Vec<f64>,
    velocities: Vec<f64>,
}

impl MoleculeDynamics {
    fn take_step(&mut self, dt: f64) -> Result<()> {
        todo!()
    }

    fn update_positions(&self, forces: &[f64], positions: &mut [f64], dt: f64) -> Result<()> {
        todo!()
    }

    fn update_velocities(&self, forces: &[f64], velocities: &mut [f64], dt: f64) -> Result<()> {
        todo!()
    }
}

fn take_md_step(
    forces: &[f64],             // F(n)
    velocities: &[f64],         // V(n)
    positions: &[f64],          // positions of atoms
    masses: &[f64],             // masses of atoms
    dt: f64,                    // Δt
    fx: F,                      // evaluate forces at new positions
    positions_new: &mut [f64],  // new positions to be updated
    velocities_new: &mut [f64], // new velocities to be updated
) where
    F: FnMut(&[f64], &mut [f64]),
{
    let m = masses.to_vector();
    let f = forces.to_vector();
    let v = velocities.to_vector();
    let r = positions.to_vector();

    // update positions
    let dr = v * dt + 0.5 * f / m * dt.powi(2);
    pot.step_toward(dr.as_slice());

    // update velecities
    let f_new = -pot.gradient().to_vector();
    let v_new = dr / dt + 0.5 * f_new / m * dt;

    todo!()
}
// e0fc00ce ends here
