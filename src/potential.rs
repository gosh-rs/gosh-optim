// [[file:../optim.note::*docs][docs:2]]
/// Represents an optimization problem with cache for avoiding unnecessary
/// function re-evaluations.
///
/// # Examples
///
/// ```ignore
/// let mut x = vec![0.0; 5];
/// let mut pot = Dynamics::new(&x, f);
/// let d = [0.2; 5];
/// pot.step_toward(&d);
/// let energy = pot.get_energy()?;
/// let forces = pot.get_forces()?;
/// let energy_old = pot.get_last_energy();
/// let forces_old = pot.get_last_forces();
/// pot.revert();
/// ```
// docs:2 ends here

// [[file:../optim.note::cc2b4eb6][cc2b4eb6]]
use super::*;
// cc2b4eb6 ends here

// [[file:../optim.note::9e96c6e5][9e96c6e5]]
#[derive(Debug, Clone)]
struct State {
    x: Vec<f64>,
    energy: Option<f64>,
    forces: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
pub struct Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>, // positions, forces => energy
{
    f: E,
    state: State,

    epsilon: f64,
    neval: usize,

    // cache previous point
    last_state: Option<State>,
}

impl State {
    fn new(x: &[f64]) -> Self {
        Self {
            x: x.to_vec(),
            energy: None,
            forces: None,
        }
    }
}
// 9e96c6e5 ends here

// [[file:../optim.note::51dd1513][51dd1513]]
impl<E> Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    /// evaluate energy and gradient at current position
    fn eval(&mut self) -> Result<(f64, &[f64])> {
        let n = self.state.x.len();

        let forces: &mut [f64] = self.state.forces.get_or_insert(vec![0.0; n]);
        let x = &self.state.x;
        let v = (self.f)(x, forces)?;
        self.state.energy = Some(v);

        self.neval += 1;

        Ok((v, forces))
    }
}
// 51dd1513 ends here

// [[file:../optim.note::1a2ff40a][1a2ff40a]]
impl<E> Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    /// Construct a Dynamics
    ///
    /// # Parameters
    ///
    /// * x: initial position
    /// * f: a closure for evaluation of energy and forces at current position.
    pub fn new(x: &[f64], f: E) -> Self {
        Dynamics {
            f,
            epsilon: 1e-8,
            neval: 0,

            state: State::new(x),
            last_state: None,
        }
    }

    /// Set epsilon for structural difference.
    pub fn set_epsilon(&mut self, eps: f64) {
        assert!(eps.is_sign_positive(), "invalid eps: {:?}", eps);
        self.epsilon = eps;
    }

    /// The number of function calls
    pub fn ncalls(&self) -> usize {
        self.neval
    }

    /// Return energy at current position if ok.
    ///
    /// The function will be evaluated when necessary.
    pub fn get_energy(&mut self) -> Result<f64> {
        match self.state.energy {
            // found cached value.
            Some(v) => Ok(v),
            // first time calculation
            None => {
                let (energy, _) = self.eval()?;
                Ok(energy)
            }
        }
    }

    /// Return energy at previous position if any
    pub fn get_last_energy(&self) -> Option<f64> {
        let state = self.last_state.as_ref()?;
        state.energy
    }

    /// Return a reference to current forces.
    ///
    /// The potential will be evaluated when necessary.
    pub fn get_forces(&mut self) -> Result<&[f64]> {
        match self.state.forces {
            // found cached value.
            Some(ref forces) => Ok(forces),
            // first time calculation
            None => {
                let (_, forces) = self.eval()?;
                Ok(forces)
            }
        }
    }

    /// Return a reference to forces at previous point.
    pub fn get_last_forces(&self) -> Option<&[f64]> {
        let state = self.last_state.as_ref()?;
        let forces = state.forces.as_ref()?;
        Some(forces)
    }

    /// Return a reference to current positions.
    pub fn positions(&self) -> &[f64] {
        &self.state.x
    }

    /// Return a reference to positions of the previous point.
    pub fn get_last_positions(&self) -> Option<&[f64]> {
        let state = self.last_state.as_ref()?;
        let x = state.x.as_ref();
        Some(x)
    }

    /// Update positions `x` with a prescribed displacement.
    ///
    /// x += displ
    pub fn step_toward(&mut self, displ: &[f64]) {
        // position changed
        if displ.vec2norm() > self.epsilon {
            // update previous point
            self.last_state = self.state.clone().into();

            // update position vector with the displacement
            self.state.x.vecadd(displ, 1.0);
            self.state.energy = None;
            self.state.forces = None;
        }
    }

    /// Revert to previous point if possible
    pub fn revert(&mut self) {
        if let Some(ref state) = self.last_state {
            self.state.clone_from(state);
        }
    }
}
// 1a2ff40a ends here

// [[file:../optim.note::aba130a2][aba130a2]]
#[test]
fn test_dynamics() -> Result<()> {
    use approx::*;

    const N: usize = 2;
    let mut x = [0.0; N];
    // f(x1, x2) = x1^2 + x2^2
    let f = |x: &[f64], f: &mut [f64]| {
        for i in 0..N {
            f[i] = -2.0 * x[i];
        }
        Ok(x.iter().map(|v| v.powi(2)).sum())
    };

    let mut pot = Dynamics::new(&x, f);
    let fx = pot.get_energy()?;
    assert_relative_eq!(fx, 0.0, epsilon = 1e-5);
    assert!(pot.get_last_energy().is_none());

    let d = [1.0, 2.0];
    pot.step_toward(&d);
    assert_eq!(pot.ncalls(), 1);
    assert_eq!(pot.get_last_energy().unwrap(), 0.0);

    let fx = pot.get_energy()?;
    assert_relative_eq!(fx, 5.0, epsilon = 1e-5);
    assert_eq!(pot.ncalls(), 2);
    let f = pot.get_forces()?;
    assert_relative_eq!(f[0], -2.0, epsilon = 1e-5);
    assert_relative_eq!(f[1], -4.0, epsilon = 1e-5);
    assert_eq!(pot.ncalls(), 2);
    assert_eq!(pot.get_last_energy().unwrap(), 0.0);
    assert_eq!(pot.get_last_forces().unwrap()[0], 0.0);

    let d = [1.0, 1.0];
    pot.step_toward(&d);
    assert_eq!(pot.positions(), &[2.0, 3.0]);

    pot.revert();
    assert_eq!(pot.positions(), &[1.0, 2.0]);
    assert_eq!(pot.get_energy()?, 5.0);
    assert_eq!(pot.ncalls(), 2);
    assert_eq!(pot.get_last_energy().unwrap(), 5.0);

    Ok(())
}
// aba130a2 ends here
