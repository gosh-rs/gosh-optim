// [[file:../optim.note::aba130a2][aba130a2]]
use gosh_core::*;
use gut::prelude::*;

use gosh_optim::{Dynamics, PotentialOutput};

#[test]
fn test_dynamics() -> Result<()> {
    use vecfx::approx::*;

    const N: usize = 2;
    let mut x = [0.0; N];
    // f(x1, x2) = x1^2 + x2^2
    let f = |x: &[f64], f: &mut [f64]| {
        for i in 0..N {
            f[i] = -2.0 * x[i];
        }
        let fx = x.iter().map(|v| v.powi(2)).sum();
        Ok(fx)
    };

    let mut pot = Dynamics::new(&x, f);
    let fx = pot.get_energy()?;
    assert_relative_eq!(fx, 0.0, epsilon = 1e-5);

    let d = [1.0, 2.0];
    pot.step_toward(&d);
    assert_eq!(pot.ncalls(), 1);

    let fx = pot.get_energy()?;
    assert_relative_eq!(fx, 5.0, epsilon = 1e-5);
    assert_eq!(pot.ncalls(), 2);
    let f = pot.get_force()?;
    assert_relative_eq!(f[0], -2.0, epsilon = 1e-5);
    assert_relative_eq!(f[1], -4.0, epsilon = 1e-5);
    assert_eq!(pot.ncalls(), 2);

    let d = [1.0, 1.0];
    pot.set_position(&d);
    assert_eq!(pot.position(), &[1.0, 1.0]);

    pot.set_epsilon(0.1);
    let d = [1.01, 1.01];
    pot.set_position(&d);
    assert_eq!(pot.position(), &[1.0, 1.0]);

    Ok(())
}
// aba130a2 ends here
