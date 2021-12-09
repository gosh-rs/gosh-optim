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
/// let fx = pot.energy();
/// let gx = pot.gradient();
/// let fx_old = pot.get_last_energy();
/// let gx_old = pot.get_last_gradient();
/// pot.revert();
/// ```
// docs:2 ends here

// [[file:../optim.note::cc2b4eb6][cc2b4eb6]]
use super::*;
// cc2b4eb6 ends here

// [[file:../optim.note::*base][base:1]]
#[derive(Debug, Clone)]
struct State {
    x: Vec<f64>,
    fx: Option<f64>,
    gx: Option<Vec<f64>>,
}

#[derive(Clone, Debug)]
pub struct Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64, // positions, forces => energy
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
            fx: None,
            gx: None,
        }
    }
}
// base:1 ends here

// [[file:../optim.note::51dd1513][51dd1513]]
impl<E> Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64,
{
    /// evaluate energy and gradient at current position
    fn eval(&mut self) -> (f64, &[f64]) {
        let n = self.state.x.len();

        let gx: &mut [f64] = self.state.gx.get_or_insert(vec![0.0; n]);
        let x = &self.state.x;
        let v = (self.f)(x, gx);
        self.state.fx = Some(v);

        self.neval += 1;

        (v, gx)
    }
}
// 51dd1513 ends here

// [[file:../optim.note::1a2ff40a][1a2ff40a]]
impl<E> Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64,
{
    /// Construct a Dynamics
    ///
    /// # Parameters
    ///
    /// * x: initial position
    /// * f: a closure for function evaluation of energy and gradient.
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

    /// Return energy at current position.
    ///
    /// The function will be evaluated when necessary.
    pub fn energy(&mut self) -> f64 {
        match self.state.fx {
            // found cached value.
            Some(v) => v,
            // first time calculation
            None => {
                let (fx, _) = self.eval();
                fx
            }
        }
    }

    /// Return energy at previous positions
    pub fn get_last_energy(&self) -> Option<f64> {
        let state = self.last_state.as_ref()?;
        state.fx
    }

    /// Return a reference to gradient at current position.
    ///
    /// The function will be evaluated when necessary.
    pub fn gradient(&mut self) -> &[f64] {
        match self.state.gx {
            // found cached value.
            Some(ref gx) => gx,
            // first time calculation
            None => {
                let (_, gx) = self.eval();
                gx
            }
        }
    }

    pub fn get_last_gradient(&self) -> Option<&[f64]> {
        let state = self.last_state.as_ref()?;
        let gx = state.gx.as_ref()?;
        Some(gx)
    }

    /// Return a reference to current position vector.
    pub fn position(&self) -> &[f64] {
        &self.state.x
    }

    /// Return a reference to last evaluated position.
    pub fn get_last_position(&self) -> Option<&[f64]> {
        let state = self.last_state.as_ref()?;
        let x = state.x.as_ref();
        Some(x)
    }

    /// Update position `x` at a prescribed displacement.
    ///
    /// x += displ
    pub fn step_toward(&mut self, displ: &[f64]) {
        // position changed
        if displ.vec2norm() > self.epsilon {
            // update previous point
            self.last_state = self.state.clone().into();

            // update position vector with the displacement
            self.state.x.vecadd(displ, 1.0);
            self.state.fx = None;
            self.state.gx = None;
        }
    }

    /// Revert to previous point if possible
    pub fn revert(&mut self) {
        if let Some(ref state) = self.last_state {
            self.state.clone_from(state);
        }

        // if let Some(ref prev) = self.x_prev {
        //     // revert position
        //     self.x.clone_from(prev);

        //     // also revert cached results
        //     self.fx = self.fx_prev;
        //     if let Some(ref prev) = self.gx_prev {
        //         if let Some(gx) = self.gx.as_mut() {
        //             gx.clone_from(prev);
        //         }
        //     }
        // }
    }
}
// 1a2ff40a ends here

// [[file:../optim.note::aba130a2][aba130a2]]
#[test]
fn test_dynamics() {
    use approx::*;

    const N: usize = 2;
    let mut x = [0.0; N];
    // f(x, y) = x^2 + y^2
    let f = |x: &[f64], g: &mut [f64]| {
        for i in 0..N {
            g[i] = 2.0 * x[i];
        }
        x.iter().map(|v| v.powi(2)).sum()
    };

    let mut pot = Dynamics::new(&x, f);
    let fx = pot.energy();
    assert_relative_eq!(fx, 0.0, epsilon = 1e-5);
    assert!(pot.get_last_energy().is_none());

    let d = [1.0, 2.0];
    pot.step_toward(&d);
    assert_eq!(pot.ncalls(), 1);
    assert_eq!(pot.get_last_energy().unwrap(), 0.0);

    let fx = pot.energy();
    assert_relative_eq!(fx, 5.0, epsilon = 1e-5);
    assert_eq!(pot.ncalls(), 2);
    let gx = pot.gradient();
    assert_relative_eq!(gx[0], 2.0, epsilon = 1e-5);
    assert_relative_eq!(gx[1], 4.0, epsilon = 1e-5);
    assert_eq!(pot.ncalls(), 2);
    assert_eq!(pot.get_last_energy().unwrap(), 0.0);
    assert_eq!(pot.get_last_gradient().unwrap()[0], 0.0);

    let d = [1.0, 1.0];
    pot.step_toward(&d);
    assert_eq!(pot.position(), &[2.0, 3.0]);

    pot.revert();
    assert_eq!(pot.position(), &[1.0, 2.0]);
    assert_eq!(pot.energy(), 5.0);
    assert_eq!(pot.ncalls(), 2);
    assert_eq!(pot.get_last_energy().unwrap(), 5.0);
}
// aba130a2 ends here

// [[file:../optim.note::709a65af][709a65af]]
#[test]
fn test_na() {
    use na::Matrix;
    use vecfx::nalgebra as na;
    let a = [1.0, 2.0];
    let b = na::DVectorSlice::from(a);
    // let c = na::DVectorSlice::from_slice(&[2.0, 3.0], a.len());
    // dbg!(b + c);
}
// 709a65af ends here
