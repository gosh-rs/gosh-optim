// [[file:../optim.note::*docs][docs:2]]
/// Represents an optimization problem with cache for avoiding unnecessary
/// function re-evaluations.
///
/// # Examples
///
/// ```ignore
/// let mut x = vec![0.0; 5];
/// let mut dyn = Dynamics::new(&x, f);
/// let d = [0.2; 5];
/// pot.step_toward(&d);
/// let fx = pot.energy();
/// let gx = pot.gradient();
/// let fx_old = pot.get_energy_prev();
/// let gx_old = pot.get_gradient_prev();
/// pot.revert();
/// ```
// docs:2 ends here

// [[file:../optim.note::cc2b4eb6][cc2b4eb6]]
use super::*;
use vecfx::*;
// cc2b4eb6 ends here

// [[file:../optim.note::1a2ff40a][1a2ff40a]]
#[derive(Clone, Debug)]
pub struct Dynamics<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64, // positions, forces => energy
{
    f: E,
    x: Vec<f64>,
    fx: Option<f64>,
    gx: Option<Vec<f64>>,

    epsilon: f64,
    neval: usize,

    // cache previous point
    x_prev: Vec<f64>,
    fx_prev: Option<f64>,
    gx_prev: Option<Vec<f64>>,
}

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
            x: x.to_vec(),
            fx: None,
            gx: None,
            neval: 0,
            x_prev: x.to_vec(),
            fx_prev: None,
            gx_prev: None,
        }
    }

    /// The number of function calls
    pub fn ncalls(&self) -> usize {
        self.neval
    }

    /// evaluate energy and gradient at current position
    fn eval(&mut self) -> (f64, &[f64]) {
        let n = self.x.len();

        let gx: &mut [f64] = self.gx.get_or_insert(vec![0.0; n]);
        let v = (self.f)(&self.x, gx);
        self.fx = Some(v);
        self.neval += 1;

        // update previous point
        self.fx_prev = self.fx;
        self.x_prev = self.x.clone();
        self.gx_prev = Some(gx.to_vec());

        (v, gx)
    }

    /// Return energy at current position.
    ///
    /// The function will be evaluated when necessary.
    pub fn energy(&mut self) -> f64 {
        match self.fx {
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
    pub fn get_energy_prev(&self) -> Option<f64> {
        self.fx_prev
    }

    /// Return a reference to gradient at current position.
    ///
    /// The function will be evaluated when necessary.
    pub fn gradient(&mut self) -> &[f64] {
        match self.gx {
            // found cached value.
            Some(ref gx) => gx,
            // first time calculation
            None => {
                let (_, gx) = self.eval();
                gx
            }
        }
    }

    pub fn get_gradient_prev(&self) -> Option<&[f64]> {
        self.gx_prev.as_ref()
    }

    /// Return a reference to current position vector.
    pub fn position(&self) -> &[f64] {
        &self.x
    }

    /// Update position `x` at a prescribed displacement.
    ///
    /// x += displ
    pub fn step_toward(&mut self, displ: &[f64]) {
        // position changed
        if displ.vec2norm() > self.epsilon {
            // update position vector with the displacement
            self.x.vecadd(displ, 1.0);
            self.fx = None;
            self.gx = None;
        }
    }

    /// Revert to previous point
    pub fn revert(&mut self) {
        self.fx = self.fx_prev;
        self.x.veccpy(&self.x_prev);
        self.gx = self.gx_prev.clone();
    }
}
// 1a2ff40a ends here

// [[file:../optim.note::aba130a2][aba130a2]]
#[test]
fn test_potential_walter() {
    const N: usize = 2;
    let mut x = [0.0; N];
    let f = |x: &[f64], g: &mut [f64]| {
        for i in 0..N {
            *g[i] = 2 * x[i];
        }
        x.iter().map(|v| v.powi(2)).sum()
    };

    let mut pot = Dynamics::new(&x, f);
    let d = [1.0, 0.0];
    pot.step_toward(&d);

    let fx = pot.energy();
    dbg!(fx);
    let gx = pot.gradient();
    dbg!(gx);

    let fx_old = pot.get_energy_prev();
    let gx_old = pot.get_gradient_prev();
    dynamics.revert();
}
// aba130a2 ends here
