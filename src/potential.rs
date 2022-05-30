// [[file:../optim.note::9f093b0e][9f093b0e]]
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
/// ```
// 9f093b0e ends here

// [[file:../optim.note::cc2b4eb6][cc2b4eb6]]
use super::*;
// cc2b4eb6 ends here

// [[file:../optim.note::0aa9588b][0aa9588b]]
#[derive(Debug, Clone)]
/// A helper struct represents the output data required for dynamics simulation.
pub struct PotentialOutput {
    /// evaluated potential energy
    pub energy: f64,
    /// evaluated force, the negative of the gradient of potential
    pub force: Vec<f64>,
}

/// Trait for potential evaluation in dynamics simulation
pub trait EvaluatePotential<U> {
    fn evaluate(&mut self, position: &[f64], output: &mut PotentialOutput) -> Result<U>;
}

impl<T> EvaluatePotential<()> for T
where
    T: FnMut(&[f64], &mut [f64]) -> Result<f64>, // position, force => energy
{
    fn evaluate(&mut self, position: &[f64], output: &mut PotentialOutput) -> Result<()> {
        let energy = self(position, &mut output.force)?;
        output.energy = energy;
        Ok(())
    }
}
// 0aa9588b ends here

// [[file:../optim.note::9e96c6e5][9e96c6e5]]
#[derive(Debug, Clone)]
struct State {
    position: Vec<f64>,
    evaluated: Option<PotentialOutput>,
}

impl State {
    fn new(x: &[f64]) -> Self {
        Self {
            position: x.to_vec(),
            evaluated: None,
        }
    }
}

/// A potential walker for dynamic simulation
pub struct Dynamics<'a, U> {
    f: Box<dyn EvaluatePotential<U> + 'a>,

    state: State,
    // cache previous point
    epsilon: f64,
    neval: usize,

    // user returned data in `evaluate` method of `EvaluatePotential` trait
    user_data: Option<U>,
}
// 9e96c6e5 ends here

// [[file:../optim.note::c39f75c1][c39f75c1]]
impl<'a, U> Dynamics<'a, U> {
    /// evaluate potential at current position
    fn eval(&mut self) -> Result<&PotentialOutput> {
        let n = self.state.position.len();
        let evaluated = self.state.evaluated.get_or_insert(PotentialOutput {
            energy: std::f64::NAN,
            force: vec![0.0; n],
        });
        let extra = self.f.evaluate(&self.state.position, evaluated)?;
        self.user_data = extra.into();
        self.neval += 1;

        Ok(evaluated)
    }

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
    pub fn new(x: &[f64], f: impl EvaluatePotential<U> + 'a) -> Self {
        Self {
            f: Box::new(f),
            epsilon: 1e-8,
            neval: 0,

            state: State::new(x),
            user_data: None,
        }
    }

    /// Return extra data returned in `evaluate` method of `EvaluatePotential`
    /// trait.
    ///
    /// The potential will be evaluated when necessary.
    pub fn get_extra(&mut self) -> Result<&U> {
        match self.user_data {
            // found cached value.
            Some(ref v) => Ok(&v),
            // first time calculation
            None => {
                let _ = self.eval()?;
                Ok(self.user_data.as_ref().unwrap())
            }
        }
    }

    /// Return energy at current position if ok.
    ///
    /// The function will be evaluated when necessary.
    pub fn get_energy(&mut self) -> Result<f64> {
        match self.state.evaluated.as_ref() {
            // found cached value.
            Some(v) => Ok(v.energy),
            // first time calculation
            None => {
                let e = self.eval()?;
                Ok(e.energy)
            }
        }
    }

    /// Return a reference to current force.
    ///
    /// The potential will be evaluated when necessary.
    pub fn get_force(&mut self) -> Result<&[f64]> {
        match self.state.evaluated {
            // found cached value.
            Some(ref e) => Ok(&e.force),
            // first time calculation
            None => {
                let e = self.eval()?;
                Ok(&e.force)
            }
        }
    }

    /// Return a reference to current position.
    pub fn position(&self) -> &[f64] {
        &self.state.position
    }

    /// The number of function calls
    pub fn ncalls(&self) -> usize {
        self.neval
    }

    /// Reset counter for potential evaluations to zero.
    pub fn recount(&mut self) {
        self.neval = 0;
    }
}
// c39f75c1 ends here

// [[file:../optim.note::1a2ff40a][1a2ff40a]]
impl<'a, U> Dynamics<'a, U> {
    /// Set epsilon for determining if structure has any substantial changes. If
    /// so, the potential will be re-evaluated automatically.
    pub fn set_epsilon(&mut self, eps: f64) {
        assert!(eps.is_sign_positive(), "invalid eps: {:?}", eps);
        self.epsilon = eps;
    }

    /// The threshold for updating current position.
    pub fn epsilon(&self) -> f64 {
        self.epsilon
    }

    /// Update position `x` with a prescribed displacement.
    ///
    /// x += displacement
    pub fn step_toward(&mut self, displacement: &[f64]) {
        // position changed
        let step_size = displacement.as_vector_slice().norm();
        assert!(step_size.is_nan(), "found invalid float numbers: {displacement:?}");
        if step_size > self.epsilon {
            // update position vector with the displacement
            self.state.position.vecadd(displacement, 1.0);
            self.state.evaluated = None;
        } else {
            info!("step size is too small: {step_size}, ignored.");
        }
    }

    /// Set current position directly.
    pub fn set_position(&mut self, position: &[f64]) {
        assert_eq!(position.len(), self.state.position.len());
        let step_size = (position.as_vector_slice() - self.state.position.as_vector_slice()).norm();
        assert!(
            step_size.is_nan(),
            "found invalid float numbers: {position:?} or {:?}",
            self.state.position
        );
        if step_size > self.epsilon {
            self.state.position.clone_from_slice(position);
            self.state.evaluated = None;
        } else {
            info!("step size is too small: {step_size}, ignored.");
        }
    }
}
// 1a2ff40a ends here

// [[file:../optim.note::b4c9a7de][b4c9a7de]]
impl<'a> Dynamics<'a, ()> {
    /// Create `Dynamics` for molecule simulation using chemical model `model`.
    pub fn from_chemical_model(
        model: &'a mut impl gosh_model::ChemicalModel,
        mut mol: gchemol::Molecule,
    ) -> Dynamics<()> {
        // handle freezing atoms/coords
        let position = mol.positions().flatten().collect_vec();
        let mask = mol.freezing_coords_mask();
        let position_opt = mask.apply(&position);
        info!("Removed {} freezing coordinates", position.len() - position_opt.len());
        Self::new(&position_opt, move |x_masked: &[f64], force: &mut [f64]| {
            let x = mask.unmask(x_masked, 0.0);
            mol.update_positions(x.as_3d().into_iter().copied());
            let mp = model.compute(&mol)?;
            let f = mp.get_forces().ok_or(format_err!("no forces"))?;
            let e = mp.get_energy().ok_or(format_err!("no energy"))?;
            // remove any contribution from freezing coords
            let f_masked = mask.apply(f.as_flat());
            assert_eq!(force.len(), f_masked.len(), "invalid mp: {mp:?}");
            force.copy_from_slice(&f_masked);

            Ok(e)
        })
    }
}
// b4c9a7de ends here
