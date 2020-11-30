// [[file:../optim.note::*imports][imports:1]]
use serde::*;

use gosh_core::gut::prelude::*;
// imports:1 ends here

// [[file:../optim.note::*base][base:1]]
#[derive(Deserialize, Debug, Clone)]
#[serde(default)]
pub(crate) struct Vars {
    pub max_step_size: f64,

    pub max_linesearch: usize,

    pub initial_step_size: f64,

    pub max_evaluations: usize,
}

impl Default for Vars {
    fn default() -> Self {
        Self {
            max_step_size: 0.1,
            initial_step_size: 1.0 / 75.0,
            max_linesearch: 1,
            max_evaluations: 0,
        }
    }
}

impl Vars {
    /// Construct `Dimer` object from environment variables
    pub fn from_env() -> Self {
        match envy::prefixed("GOSH_OPTIM_").from_env::<Self>() {
            Ok(vars) => {
                debug!("Found gosh-optim env variables.");
                vars
            }
            Err(error) => {
                warn!("No relevant environment variables found.");
                Self::default()
            }
        }
    }
}
// base:1 ends here
