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
/// let force = pot.get_force()?;
/// let energy_old = pot.get_last_energy();
/// let force_old = pot.get_last_force();
/// pot.revert();
/// ```
// docs:2 ends here

// [[file:../optim.note::cc2b4eb6][cc2b4eb6]]
use super::*;
// cc2b4eb6 ends here

// [[file:../optim.note::9e96c6e5][9e96c6e5]]
/// Trait for potential evaluation in dynamics simulation
pub trait EvaluateEnergyForce {
    /// Return energy and force of potential at current position.
    fn evaluate(&mut self, position: &[f64], force: &mut [f64]) -> Result<f64>;
}

impl<T> EvaluateEnergyForce for T
where
    T: FnMut(&[f64], &mut [f64]) -> Result<f64>, // position, force => energy
{
    fn evaluate(&mut self, position: &[f64], force: &mut [f64]) -> Result<f64> {
        self(position, force)
    }
}

#[derive(Debug, Clone)]
struct State {
    x: Vec<f64>,
    energy: Option<f64>,
    force: Option<Vec<f64>>,
}

/// A potential walker for dynamic simulation
pub struct Dynamics<'a> {
    f: Box<dyn EvaluateEnergyForce + 'a>,
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
            force: None,
        }
    }
}
// 9e96c6e5 ends here

// [[file:../optim.note::51dd1513][51dd1513]]
impl<'a> Dynamics<'a> {
    /// evaluate energy and gradient at current position
    fn eval(&mut self) -> Result<(f64, &[f64])> {
        let n = self.state.x.len();

        let force: &mut [f64] = self.state.force.get_or_insert(vec![0.0; n]);
        let x = &self.state.x;
        let v = self.f.evaluate(x, force)?;
        self.state.energy = Some(v);

        self.neval += 1;

        Ok((v, force))
    }
}
// 51dd1513 ends here

// [[file:../optim.note::1a2ff40a][1a2ff40a]]
impl<'a> Dynamics<'a> {
    /// Construct a Dynamics
    ///
    /// # Parameters
    ///
    /// * x: initial position
    /// * f: a closure for evaluation of energy and force at current position.
    ///   in closure f:
    ///      - the first paramater is the position
    ///      - the second is the force to be updated
    ///      - the return value is the energy
    pub fn new(x: &[f64], f: impl EvaluateEnergyForce + 'a) -> Self {
        Dynamics {
            f: Box::new(f),
            epsilon: 1e-8,
            neval: 0,

            state: State::new(x),
            last_state: None,
        }
    }

    /// Set epsilon for determining if structure has any substantial changes. If
    /// so, the potential will be re-evaluated automatically.
    pub fn set_epsilon(&mut self, eps: f64) {
        assert!(eps.is_sign_positive(), "invalid eps: {:?}", eps);
        self.epsilon = eps;
    }

    /// The number of function calls
    pub fn ncalls(&self) -> usize {
        self.neval
    }

    /// Reset counter for potential evaluations to zero.
    pub fn recount(&mut self) {
        self.neval = 0;
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

    /// Return a reference to current force.
    ///
    /// The potential will be evaluated when necessary.
    pub fn get_force(&mut self) -> Result<&[f64]> {
        match self.state.force {
            // found cached value.
            Some(ref force) => Ok(force),
            // first time calculation
            None => {
                let (_, force) = self.eval()?;
                Ok(force)
            }
        }
    }

    /// Return a reference to force at previous point.
    pub fn get_last_force(&self) -> Option<&[f64]> {
        let state = self.last_state.as_ref()?;
        let force = state.force.as_ref()?;
        Some(force)
    }

    /// Return a reference to current position.
    pub fn position(&self) -> &[f64] {
        &self.state.x
    }

    /// Return a reference to position of the previous point.
    pub fn get_last_position(&self) -> Option<&[f64]> {
        let state = self.last_state.as_ref()?;
        let x = state.x.as_ref();
        Some(x)
    }

    /// Update position `x` with a prescribed displacement.
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
            self.state.force = None;
        }
    }

    /// Set current position directly.
    pub fn set_position(&mut self, position: &[f64]) {
        assert_eq!(position.len(), self.state.x.len());
        self.last_state = self.state.clone().into();
        self.state.x.clone_from_slice(position);
        self.state.energy = None;
        self.state.force = None;
    }

    /// Revert to previous point if possible
    pub fn revert(&mut self) {
        if let Some(ref state) = self.last_state {
            self.state.clone_from(state);
        }
    }
}
// 1a2ff40a ends here
